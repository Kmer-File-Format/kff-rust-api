//! Declaration of trait to read and write kmer section

/* crate use */
use bitvec::prelude::*;
use byteorder::*;

/* local use */
use crate::kmer::Kmer;
use crate::seq2bits::Seq2Slice;
use crate::utils::{BitBox, BitOrd, BitVec};
use crate::*;

use crate::error::LocalResult;

/// Reader trait must be implement by Struct want read kmer KFF section
pub trait Reader<'input, R>: Iterator<Item = crate::Result<Kmer>>
where
    R: std::io::Read + 'input,
    R: ?Sized,
{
    /// Get size of kmers
    fn k(&self) -> u64;

    /// Get a reference to the read stream
    fn input(&mut self) -> &mut R;

    /// Get the max number of kmer in a block
    fn max_kmer(&self) -> u64;

    /// Get the actual number of kmer remains in block
    fn block_n(&self) -> u64;

    /// Get the size of data attach to each kmer in bytes
    fn data_size(&self) -> u64;

    /// Read on kmer block and return number of bytes read
    fn read_block(&mut self) -> crate::Result<usize>;

    /// Get the sequence of actual block
    fn block_seq(&self) -> &Seq2Slice;

    /// Get the data of actutal block
    fn block_data(&self) -> &[u8];

    /// Get bit used to perform reverse encoding
    fn rev_encoding(&self) -> u8;

    /// Reduce block_n by one
    fn decrease_n(&mut self);

    /// Read the number of kmer in block
    fn read_n(&mut self) -> crate::Result<u64> {
        if self.max_kmer() == 1 {
            Ok(1)
        } else {
            let max_value = self.max_kmer();
            utils::read_dynamic_size_field(self.input(), max_value)
        }
    }

    /// Read sequence of actual block need number of nucleotide want be to read
    fn read_seq(&mut self, nb_nuc: u64) -> crate::Result<BitBox> {
        let buf_len = utils::ceil_to_8(nb_nuc * 2) as usize / 8;
        let mut buffer = vec![0u8; buf_len];

        self.input().read_exact(&mut buffer).map_local()?;

        let bit_buffer = BitVec::from_vec(buffer);

        Ok(BitBox::from_bitslice(
            &bit_buffer[(buf_len * 8 - nb_nuc as usize * 2)..],
        ))
    }

    /// Read data of actual block
    fn read_data(&mut self) -> crate::Result<Vec<u8>> {
        let mut buffer = vec![0u8; (self.block_n() * self.data_size()) as usize];

        self.input().read_exact(&mut buffer).map_local()?;
        Ok(buffer)
    }

    /// Get the next kmer
    fn next_kmer(&mut self) -> crate::Result<Kmer> {
        if self.block_n() == 0 && self.read_block()? == 0 {
            return Ok(Kmer::new(
                bitbox![BitOrd, u8; 0; 0],
                Vec::new().into_boxed_slice(),
            ));
        }

        let index = self.block_seq().len() / 2 - (self.k() + self.block_n() - 1) as usize;

        let kmer = BitBox::from_bitslice(
            &self.block_seq()[(index * 2)..((index + self.k() as usize) * 2)],
        );

        let data = Vec::from(
            &self.block_data()
                [(index * self.data_size() as usize)..((index + 1) * self.data_size() as usize)],
        );

        self.decrease_n();

        Ok(Kmer::new(kmer, data.into_boxed_slice()))
    }
}

/// Write trait must be implement by Struct want write kmer KFF section
pub trait Writer<'output, W>: Sized
where
    W: std::io::Write + std::io::Seek + 'output,
{
    /// Get the size of data attach to each kmer in bytes
    fn data_size(&self) -> u64;

    /// Get the max number of kmer in a block
    fn max(&self) -> u64;

    /// Position of nb_block_offset need to write the number of block in section after finish to write section
    fn nb_block_offset(&self) -> u64;

    /// Number of block in this section
    fn nb_block(&self) -> u32;

    /// Return True if section is close don't write any other block in this section
    fn is_close(&self) -> bool;

    /// Get a reference of output
    fn output(&mut self) -> &mut W;

    // Setter
    /// Change number of block in section
    fn set_nb_block(&mut self, value: u32);

    /// Change state of [Writer::is_close]
    fn set_close(&mut self, value: bool);

    // Computation
    /// Compute the number of kmer from number of sequence
    fn nb_kmer(&self, seq_len: usize) -> usize;

    // Default implementation
    /// Close the section, can return an error because number of block can only by write when number of block is know
    fn close(&mut self) -> crate::Result<usize> {
        if self.is_close() {
            Ok(0)
        } else {
            let offset = self.nb_block_offset();
            let nb_block = self.nb_block();

            self.output()
                .seek(std::io::SeekFrom::Start(offset))
                .map_local()?;
            self.output()
                .write_u32::<LittleEndian>(nb_block)
                .map_local()?;
            self.output().seek(std::io::SeekFrom::End(0)).map_local()?;
            self.set_close(true);

            Ok(1)
        }
    }

    /// Increase the number of block
    fn increment_nb_block(&mut self) -> crate::Result<()> {
        match self.nb_block().checked_add(1) {
            Some(val) => {
                self.set_nb_block(val);
                Ok(())
            }
            None => Err(error::Error::Data(error::Data::ToManyBlock)),
        }
    }

    /// Verify the number of the number of kmer and data length match and not upper than max number of kmer
    fn check_block(&mut self, seq_len: usize, data_len: usize) -> crate::Result<usize> {
        let nb_kmer = self.nb_kmer(seq_len);
        let nb_data = data_len / self.data_size() as usize;

        if nb_data != nb_kmer {
            return Err(error::Error::Data(error::Data::NbKmerNbDataDiff));
        }

        if nb_kmer > self.max() as usize {
            return Err(error::Error::Data(error::Data::NUpperThanMax));
        }

        Ok(nb_kmer)
    }
}
