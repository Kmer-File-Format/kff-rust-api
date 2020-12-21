/* crate use */
use bitvec::prelude::*;
use byteorder::*;

/* local use */
use crate::kmer::Kmer;
use crate::utils::{BitBox, BitOrd, BitSlice, BitVec};
use crate::*;

pub trait Reader<'input, R>: Iterator<Item = crate::Result<Kmer>>
where
    R: std::io::Read + 'input,
    R: ?Sized,
{
    fn k(&self) -> u64;
    fn input(&mut self) -> &mut R;
    fn max_kmer(&self) -> u64;
    fn block_n(&self) -> u64;
    fn data_size(&self) -> u64;
    fn read_block(&mut self) -> crate::Result<usize>;
    fn block_seq(&self) -> &BitSlice;
    fn block_data(&self) -> &[u8];
    fn rev_encoding(&self) -> u8;
    fn decrease_n(&mut self);

    fn read_n(&mut self) -> crate::Result<u64> {
        if self.max_kmer() == 1 {
            Ok(1)
        } else {
            let max_value = self.max_kmer();
            utils::read_dynamic_size_field(self.input(), max_value)
        }
    }

    fn read_seq(&mut self, nb_nuc: u64) -> crate::Result<BitBox> {
        let buf_len = utils::ceil_to_8(nb_nuc * 2) as usize / 8;
        let mut buffer = vec![0u8; buf_len];

        self.input().read_exact(&mut buffer)?;

        let bit_buffer = BitVec::from_vec(buffer);

        Ok(BitBox::from_bitslice(
            &bit_buffer[(buf_len * 8 - nb_nuc as usize * 2)..],
        ))
    }

    fn read_data(&mut self) -> crate::Result<Vec<u8>> {
        let mut buffer = vec![0u8; (self.block_n() * self.data_size()) as usize];

        self.input().read_exact(&mut buffer)?;
        Ok(buffer)
    }

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

pub trait Writer<'output, W>: Drop + Sized
where
    W: std::io::Write + std::io::Seek + 'output,
{
    // Getter
    fn data_size(&self) -> u64;
    fn max(&self) -> u64;
    fn nb_block_offset(&self) -> u64;
    fn nb_block(&self) -> u32;

    fn is_close(&self) -> bool;

    fn output(&mut self) -> &mut W;

    // Setter
    fn set_nb_block(&mut self, value: u32);
    fn set_close(&mut self, value: bool);

    // Computation
    fn nb_kmer(&self, seq_len: usize) -> usize;

    // Default implementation
    fn close(&mut self) -> crate::Result<usize> {
        if self.is_close() {
            Ok(0)
        } else {
            let offset = self.nb_block_offset();
            let nb_block = self.nb_block();

            self.output().seek(std::io::SeekFrom::Start(offset))?;
            self.output().write_u32::<LittleEndian>(nb_block)?;
            self.output().seek(std::io::SeekFrom::End(0))?;
            self.set_close(true);

            Ok(1)
        }
    }

    fn increment_nb_block(&mut self) -> crate::Result<()> {
        match self.nb_block().checked_add(1) {
            Some(val) => {
                self.set_nb_block(val);
                Ok(())
            }
            None => Err(Box::new(error::Data::ToManyBlock)),
        }
    }

    fn check_block(&mut self, seq_len: usize, data_len: usize) -> crate::Result<usize> {
        let nb_kmer = self.nb_kmer(seq_len);
        let nb_data = data_len / self.data_size() as usize;

        if nb_data != nb_kmer {
            return Err(Box::new(error::Data::NbKmerNbDataDiff));
        }

        if nb_kmer > self.max() as usize {
            return Err(Box::new(error::Data::NUpperThanMax));
        }

        Ok(nb_kmer)
    }
}
