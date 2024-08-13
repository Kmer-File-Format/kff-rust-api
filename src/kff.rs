//! Kmer File Format Rust parser

/* std use */

/* crate use */
use std::io::Seek;

/* project use */
use crate::error;
use crate::section;
use crate::utils;
use crate::GlobalIndex;
use crate::KffRead;
use crate::Kmer;
use crate::KmerIterator;

use crate::section::values::AbcValues as _;

/// Struct to read a kff file
#[derive(getset::Getters, getset::Setters, getset::MutGetters)]
#[getset(get = "pub")]
pub struct Kff<T> {
    /// Inner read source
    inner: T,

    /// Header extract from `inner`
    #[getset(set = "pub", get_mut = "pub")]
    header: section::Header,

    /// Current Values extract from `inner`
    #[getset(set = "pub", get_mut = "pub")]
    values: section::Values,

    /// GlobalIndex present only if inner is seekable and first section is index or footer contains first_index
    index: Option<utils::GlobalIndex>,
}

impl<R> Kff<R>
where
    R: std::io::Read + std::io::BufRead + crate::KffRead,
{
    /// Create a new Kff reader by accept mutable reference on [std::io::Read]
    pub fn read(mut inner: R) -> error::Result<Self> {
        let header = section::Header::read(&mut inner)?;
        let values = section::Values::default();

        Ok(Self {
            inner,
            header,
            values,
            index: None,
        })
    }

    /// Consume Kff object to create a KmerIterator
    pub fn kmers(self) -> KmerIterator<R> {
        KmerIterator::new(self)
    }

    /// Read Kff until last kmer section
    pub fn next_kmer_section(&mut self) -> std::option::Option<error::Result<Vec<Kmer>>> {
        loop {
            match self.inner.read_u8() {
                Ok(b'v') => {
                    self.values = {
                        match section::Values::read(&mut self.inner) {
                            Err(e) => return Some(Err(e)),
                            Ok(v) => v,
                        }
                    }
                }
                Ok(b'r') => match section::Raw::new(&self.values) {
                    Ok(section) => return Some(section.read(&mut self.inner)),
                    Err(e) => return Some(Err(e)),
                },
                Ok(b'm') => match section::Minimizer::new(&self.values) {
                    Ok(section) => return Some(section.read(&mut self.inner)),
                    Err(e) => return Some(Err(e)),
                },
                Ok(b'K') => return None, // It's the begin of last signature stop reading
                Ok(b'i') => match section::Index::skip(&mut self.inner) {
                    Err(e) => return Some(Err(e)),
                    Ok(_) => continue,
                },
                Ok(e) => return Some(Err(error::Kff::NotASectionPrefix(e).into())), // Any other value is an error
                Err(e) => return Some(Err(e)),
            }
        }
    }
}

impl Kff<std::io::BufReader<std::fs::File>> {
    /// Create a new Kff by read file match with path
    pub fn open<P>(path: P) -> error::Result<Self>
    where
        P: std::convert::AsRef<std::path::Path>,
    {
        std::fs::File::open(&path)
            .map(std::io::BufReader::new)
            .map(Kff::read)?
    }

    /// Create a Kff and generate a global index
    pub fn with_index<P>(path: P) -> error::Result<Self>
    where
        P: std::convert::AsRef<std::path::Path>,
    {
        let mut inner = std::fs::File::open(&path).map(std::io::BufReader::new)?;

        let header = section::Header::read(&mut inner)?;
        let values = section::Values::default();

        let pos_first_section = inner.stream_position()?;
        let index = match utils::GlobalIndex::new(&mut inner, pos_first_section) {
            Ok(index) => Some(index),
            Err(error::Error::Kff(error::Kff::NotAnIndex)) => {
                let value = Kff::load_footer(&mut inner)?;

                Some(GlobalIndex::new(
                    &mut inner,
                    *value.get("first_index").ok_or(error::Kff::NoFirstIndex)?,
                )?)
            }
            Err(e) => return Err(e),
        };

        Ok(Self {
            inner,
            header,
            values,
            index,
        })
    }
}

impl<R> Kff<R>
where
    R: std::io::Read + std::io::Seek + KffRead,
{
    /// Check readable match with a KFF file
    pub fn check(&mut self) -> error::Result<bool> {
        let cursor_position = self.inner.stream_position()?;

        self.inner.seek(std::io::SeekFrom::Start(0))?;
        let magic_number = self.inner.read_n_bytes::<3>()?;
        if &magic_number != b"KFF" {
            return Err(error::Kff::MissingMagic("start".to_string()).into());
        }

        self.inner.seek(std::io::SeekFrom::End(-3))?;
        let magic_number = self.inner.read_n_bytes::<3>()?;
        if &magic_number != b"KFF" {
            return Err(error::Kff::MissingMagic("end".to_string()).into());
        }

        self.inner.seek(std::io::SeekFrom::Start(cursor_position))?;

        Ok(true)
    }

    /// Load footer, assume last section is a value and last value of this section is footer_size
    fn load_footer(inner: &mut R) -> error::Result<section::Values> {
        inner.seek(std::io::SeekFrom::End(-11))?;
        let footer_size = inner.read_u64()?;

        inner.seek(std::io::SeekFrom::End(-(footer_size as i64 + 3)))?;

        let v = inner.read_u8()?;
        if v != b'v' {
            Err(error::Kff::FooterSizeNotCorrect.into())
        } else {
            section::Values::read(inner)
        }
    }

    /// Get kmer of nth section in index.
    ///
    /// If index isn't set return an Error
    /// If we didn't found section value before target section return an Error
    /// If section isn't a kmer section return an Error
    pub fn kmer_of_section(&mut self, n: usize) -> error::Result<Vec<Kmer>> {
        let index = self
            .index
            .as_ref()
            .ok_or(error::Error::Kff(error::Kff::NoIndex))?;

        self.values = match index.pair()[..n].iter().rev().find(|x| x.0 == b'v') {
            Some((_t, p)) => {
                self.inner.seek(std::io::SeekFrom::Start(p + 1))?;
                section::Values::read(&mut self.inner)?
            }
            None => return Err(error::Kff::NoValueSectionBeforeTarget.into()),
        };

        self.inner
            .seek(std::io::SeekFrom::Start(index.pair()[n].1))?;
        match self.inner.read_u8()? {
            b'r' => match section::Raw::new(&self.values) {
                Ok(section) => section.read(&mut self.inner),
                Err(e) => Err(e),
            },
            b'm' => match section::Minimizer::new(&self.values) {
                Ok(section) => section.read(&mut self.inner),
                Err(e) => Err(e),
            },
            _ => Err(error::Kff::NotAKmerSection.into()),
        }
    }
}

impl<T> std::io::Seek for Kff<T>
where
    T: std::io::Read + std::io::Seek,
{
    fn seek(&mut self, pos: std::io::SeekFrom) -> std::io::Result<u64> {
        self.inner.seek(pos)
    }
}

impl<W> Kff<W>
where
    W: std::io::Write + crate::KffWrite,
{
    /// Create a Kff object to write in inner
    pub fn write(mut inner: W, header: section::Header) -> error::Result<Self> {
        header.write(&mut inner)?;

        Ok(Self {
            inner,
            header,
            values: section::Values::default(),
            index: None,
        })
    }

    /// Write a Values section
    pub fn write_values(&mut self, values: section::Values) -> error::Result<()> {
        self.values = values;

        self.inner.write_bytes(b"v")?;
        self.values.write(&mut self.inner)?;

        Ok(())
    }

    /// Write a Index section
    pub fn write_index(&mut self, index: section::Index) -> error::Result<()> {
        self.inner.write_bytes(b"i")?;
        index.write(&mut self.inner)
    }

    /// Write a Raw section
    pub fn write_raw(
        &mut self,
        raw: section::Raw,
        blocks: &[section::block::Block],
    ) -> error::Result<()> {
        self.inner.write_bytes(b"r")?;
        raw.write(&mut self.inner, blocks)
    }

    /// Write a Minimizer section
    pub fn write_minimizer(
        &mut self,
        section: section::Minimizer,
        minimizer: crate::Seq2Bit,
        blocks: &[section::block::Block],
    ) -> error::Result<()> {
        self.inner.write_bytes(b"m")?;
        section.write(&mut self.inner, minimizer, blocks)
    }

    /// Finalize write the final signature
    pub fn finalize(&mut self) -> error::Result<()> {
        self.inner.write_bytes(b"KFF")?;
        self.inner.flush()?;

        Ok(())
    }
}

impl Kff<std::io::BufWriter<std::fs::File>> {
    /// Intialize a Kff object to write file
    pub fn create<P>(path: P, header: section::Header) -> error::Result<Self>
    where
        P: std::convert::AsRef<std::path::Path>,
    {
        std::fs::File::create(&path)
            .map(std::io::BufWriter::new)
            .map(|x| Kff::write(x, header))?
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Read;
    use std::io::Seek;

    const KFF_FILE: &[u8] = &[
        b'K', b'F', b'F',       // Magic number
        1,          // Major version number
        0,          // Minor version number
        0b00101110, // Encoding
        1,          // Uniq kmer
        0,          // Canonical kmer
        0, 0, 0, 4, // Free space size
        b't', b'e', b's', b't', // Free space
        b'v', // Footer
        0, 0, 0, 0, 0, 0, 0, 1, // Footer nb variables
        b'f', b'o', b'o', b't', b'e', b'r', b'_', b's', b'i', b'z', b'e',
        0, // name of variable footer_size
        0, 0, 0, 0, 0, 0, 0, 29, // value of variable footer_size
        b'K', b'F', b'F', // Magic number
    ];

    use std::io::Write;

    #[test]
    fn create_kff_reader() -> error::Result<()> {
        let inner = std::io::Cursor::new(KFF_FILE.to_vec());

        let mut reader = Kff::read(inner)?;
        assert!(reader.check()?);

        let mut tmpfile = tempfile::NamedTempFile::new()?;
        tmpfile.write_all(KFF_FILE)?;
        let mut reader: Kff<std::io::BufReader<std::fs::File>> =
            Kff::<std::io::BufReader<std::fs::File>>::open(tmpfile.path())?;

        assert!(reader.check()?);

        Ok(())
    }

    #[test]
    fn check_header() -> error::Result<()> {
        let inner = std::io::Cursor::new(KFF_FILE.to_vec());

        let reader = Kff::read(inner)?;

        assert_eq!(reader.header().major_version(), &1);
        assert_eq!(reader.header().minor_version(), &0);
        assert_eq!(reader.header().encoding(), &0b00101110);
        assert_eq!(reader.header().uniq_kmer(), &true);
        assert_eq!(reader.header().canonical_kmer(), &false);
        assert_eq!(reader.header().free_block(), b"test");

        Ok(())
    }

    #[test]
    fn check() -> error::Result<()> {
        let inner_len = KFF_FILE.len();
        let mut readable = std::io::Cursor::new(KFF_FILE.to_vec());

        let mut file = Kff::read(readable.clone())?;

        assert!(file.check()?); // Header init and check work

        // move cursor before check
        let mut second = readable.clone();
        let mut file = Kff::read(second.clone())?;
        second.seek_relative(10)?;
        assert!(matches!(second.stream_position(), Ok(10)));
        assert!(file.check()?);
        assert!(matches!(second.stream_position(), Ok(10)));

        readable.get_mut()[1] = b'K';
        let file = Kff::read(readable.clone());
        assert!(file.is_err()); // Header init failled

        readable.get_mut()[1] = b'F';
        readable.get_mut()[inner_len - 1] = b'K';
        let mut file = Kff::read(readable.clone())?;

        assert!(file.check().is_err()); // Header init work but footer failled

        Ok(())
    }

    #[test]
    fn load_footer() -> error::Result<()> {
        let inner_len = KFF_FILE.len();
        let mut inner = std::io::Cursor::new(KFF_FILE.to_vec());

        let mut truth = section::Values::new();
        truth.insert("footer_size".to_string(), 29);

        assert_eq!(Kff::load_footer(&mut inner)?, truth);

        inner.get_mut()[inner_len - 32] = b'f';

        assert!(Kff::load_footer(&mut inner).is_err());

        Ok(())
    }

    #[test]
    fn seek() -> error::Result<()> {
        let inner = std::io::Cursor::new(KFF_FILE.to_vec());
        let mut reader = Kff::read(inner)?;

        assert_eq!(reader.seek(std::io::SeekFrom::Start(2))?, 2);

        Ok(())
    }

    #[test]
    fn write() -> error::Result<()> {
        let file = tempfile::NamedTempFile::new()?;
        let header = section::Header::new(1, 0, 0b00011011, true, true, b"".to_vec())?;
        let mut writer = Kff::create(file.path(), header)?;

        let mut values = section::Values::default();
        values.insert("k".to_string(), 5);
        values.insert("m".to_string(), 3);
        values.insert("ordered".to_string(), false as u64);
        values.insert("max".to_string(), 200);
        values.insert("data_size".to_string(), 1);

        writer.write_values(values.clone())?;

        writer.write_raw(section::Raw::new(&values)?, &[
	    section::block::Block {
                k: 5,
                data_size: 1,
                kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
				vec![1, 2, 3]),
		minimizer_offset: 0,
		offset: 0,
            },
            section::block::Block{
                k: 5,
                data_size: 1,
                kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
				vec![1, 2]),
		minimizer_offset: 0,
		offset: 0,
            },
            section::block::Block {
                k: 5,
		data_size: 1,
                kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
				vec![1]),
		minimizer_offset: 0,
		offset: 0,
            },
	])?;

        writer.write_minimizer(
	    section::Minimizer::new(&values)?,
	    bitvec::bitbox![u8, bitvec::order::Msb0; 0, 1, 1, 0, 1, 1],
            &[
                section::block::Block{
                    k: 5,
                    data_size: 1,
                    kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                     vec![1, 2, 3]),
		    minimizer_offset: 1,
		    offset: 0,
                },
                section::block::Block {
		    k: 5,
		    data_size: 1,
		    kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1],
                    vec![1, 2]),
		    minimizer_offset: 1,
		    offset: 0,
                },
                section::block::Block {
		    k: 5,
                    data_size: 1,
		    kmer: Kmer::new(bitvec::bitbox![u8, bitvec::order::Msb0; 0, 0, 0, 1, 1, 0, 1, 1, 1, 1],
		    vec![1]),
		    minimizer_offset: 1,
		    offset: 0,
                }
            ],
	)?;

        writer.write_index(section::Index::new(
            vec![(b'v', -30), (b'r', -25), (b'm', -20)],
            0,
        ))?;

        writer.finalize()?;

        let mut inner = Vec::new();
        let (_, path) = file.keep().unwrap();
        let mut t = std::fs::File::open(path)?;
        t.read_to_end(&mut inner)?;

        assert_eq!(
            inner,
            vec![
                b'K', b'F', b'F', //
                1, 0,  // Version number
                27, // Encoding
                1, 1, // Uniq, Canonical
                0, 0, 0, 0, // Free space size length
                b'v', 0, 0, 0, 0, 0, 0, 0, 5, // Five values
                b'o', b'r', b'd', b'e', b'r', b'e', b'd', 0, 0, 0, 0, 0, 0, 0, 0, 0, //
                b'd', b'a', b't', b'a', b'_', b's', b'i', b'z', b'e', 0, 0, 0, 0, 0, 0, 0, 0, 1,
                b'm', 0, 0, 0, 0, 0, 0, 0, 0, 3, //
                b'k', 0, 0, 0, 0, 0, 0, 0, 0, 5, //
                b'm', b'a', b'x', 0, 0, 0, 0, 0, 0, 0, 0, 200, //
                b'r', 0, 0, 0, 0, 0, 0, 0, 3, // Three block
                3, 27, 244, 1, 2, 3, // Three kmer in block
                2, 27, 240, 1, 2, // Two kmer in block
                1, 27, 192, 1,    // One kmer in block
                b'm', //
                108,  // minimizer sequence
                0, 0, 0, 0, 0, 0, 0, 3, // Three block
                3, 1, 61, 1, 2, 3, // Three kmer minimizer at offset 1
                2, 1, 60, 1, 2, // Two kmer minimizer at offset 1
                1, 1, 48, 1, // One kmer minimizer at offset 1
                b'i', 0, 0, 0, 0, 0, 0, 0, 3, // Three section indexed
                b'v', 255, 255, 255, 255, 255, 255, 255, 226, // Value section
                b'r', 255, 255, 255, 255, 255, 255, 255, 231, // Raw section
                b'm', 255, 255, 255, 255, 255, 255, 255, 236, // Minimizer section
                0, 0, 0, 0, 0, 0, 0, 0, // No other index
                b'K', b'F', b'F', //
            ]
        );

        Ok(())
    }
}
