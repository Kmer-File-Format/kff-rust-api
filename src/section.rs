//! Types shared by the reader and the writer: section types, variables and index sections

/* std use */
use std::collections::BTreeMap;

/// Variables of a values section ('v'), sorted by name.
///
/// The C++ API stores them in a `std::map`: the variables of a section are written in the order of
/// their names, and this crate does the same.
pub type Values = BTreeMap<String, u64>;

/// Type of a section, given by its first byte
#[derive(Clone, Copy, std::fmt::Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum SectionType {
    /// Values declarations, 'v'
    Values,
    /// Index, 'i'
    Index,
    /// Raw sequences, 'r'
    Raw,
    /// Sequences with minimizers, 'm'
    Minimizer,
    /// Any other byte
    Other(u8),
}

impl SectionType {
    /// Section type of the first byte of a section
    pub fn from_byte(byte: u8) -> Self {
        match byte {
            b'v' => Self::Values,
            b'i' => Self::Index,
            b'r' => Self::Raw,
            b'm' => Self::Minimizer,
            other => Self::Other(other),
        }
    }

    /// First byte of a section of this type
    pub fn as_byte(self) -> u8 {
        match self {
            Self::Values => b'v',
            Self::Index => b'i',
            Self::Raw => b'r',
            Self::Minimizer => b'm',
            Self::Other(other) => other,
        }
    }

    /// True for the sections that contain kmers (raw and minimizer sections)
    pub fn has_kmers(self) -> bool {
        matches!(self, Self::Raw | Self::Minimizer)
    }
}

impl std::convert::From<u8> for SectionType {
    fn from(byte: u8) -> Self {
        Self::from_byte(byte)
    }
}

impl std::convert::From<SectionType> for u8 {
    fn from(section_type: SectionType) -> Self {
        section_type.as_byte()
    }
}

/// Content of an index section ('i')
#[derive(Clone, std::fmt::Debug, PartialEq, Eq, Default)]
pub struct IndexSection {
    /// Position of the section in the file (its 'i' byte)
    pub beginning: u64,

    /// Indexed sections: position relative to the end of the index section, and type
    pub entries: BTreeMap<i64, SectionType>,

    /// Position of the next index section relative to the end of this one, 0 for the last one
    pub next_index: i64,
}

impl IndexSection {
    /// Size of the section in bytes
    pub fn size(&self) -> u64 {
        index_section_size(self.entries.len() as u64)
    }

    /// Position of the first byte after the section, the origin of the relative positions
    pub fn end(&self) -> u64 {
        self.beginning + self.size()
    }

    /// Absolute positions in the file of the indexed sections, with their type
    pub fn absolute_entries(&self) -> impl Iterator<Item = (u64, SectionType)> + '_ {
        let end = self.end() as i64;
        self.entries
            .iter()
            .map(move |(relative, section_type)| ((end + relative) as u64, *section_type))
    }
}

/// Size of an index section with `nb_entries` entries
pub(crate) fn index_section_size(nb_entries: u64) -> u64 {
    1 + 8 + 9 * nb_entries + 8
}

/// Value of a big endian unsigned integer
#[inline]
pub(crate) fn load_big_endian(bytes: &[u8]) -> u64 {
    bytes
        .iter()
        .fold(0, |value, byte| (value << 8) | *byte as u64)
}

/// Big endian representation of value on `nb_bytes` bytes (at most 8), the highest bytes are lost
/// as in the C++ API.
#[inline]
pub(crate) fn store_big_endian(value: u64, nb_bytes: usize) -> [u8; 8] {
    let mut bytes = [0; 8];
    for (i, byte) in bytes[..nb_bytes].iter_mut().rev().enumerate() {
        *byte = (value >> (8 * i)) as u8;
    }
    bytes
}

/// Number of bytes of a field that can contain values up to `max`: `ceil(log2(max))` bits.
///
/// The C++ API computes this width with floating point numbers; the same computation is used here to
/// read and write the same bytes for any value.
#[inline]
pub(crate) fn bytes_for_max(max: u64) -> usize {
    let bits = (max as f64).log2().ceil();
    let bits = if bits > 0.0 { bits as u64 } else { 0 };
    ((bits + 7) / 8) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn section_types() {
        for byte in 0..=255u8 {
            assert_eq!(SectionType::from_byte(byte).as_byte(), byte);
        }
        assert_eq!(SectionType::from(b'm'), SectionType::Minimizer);
        assert!(SectionType::Raw.has_kmers());
        assert!(!SectionType::Values.has_kmers());
    }

    #[test]
    fn widths() {
        // Same widths as the tests of the C++ API
        let widths: &[(u64, usize)] = &[
            (1, 0),
            (2, 1),
            (3, 1),
            (255, 1),
            (256, 1),
            (257, 2),
            (65535, 2),
            (65536, 2),
            (65537, 3),
            (1 << 24, 3),
            ((1 << 24) + 1, 4),
            (1 << 32, 4),
            ((1 << 32) + 1, 5),
            (u64::MAX, 8),
        ];
        for (max, bytes) in widths {
            assert_eq!(bytes_for_max(*max), *bytes, "max {max}");
        }
    }

    #[test]
    fn big_endian() {
        assert_eq!(&store_big_endian(0x0102, 2)[..2], &[1, 2]);
        assert_eq!(&store_big_endian(0x010203, 2)[..2], &[2, 3]);
        assert_eq!(&store_big_endian(u64::MAX, 8), &[255; 8]);
        assert_eq!(load_big_endian(&[1, 2, 3]), 0x010203);
        assert_eq!(load_big_endian(&[]), 0);
    }

    #[test]
    fn index_positions() {
        let mut index = IndexSection {
            beginning: 100,
            next_index: 0,
            ..Default::default()
        };
        index.entries.insert(-50, SectionType::Values);
        index.entries.insert(10, SectionType::Raw);
        assert_eq!(index.size(), 1 + 8 + 18 + 8);
        assert_eq!(
            index.absolute_entries().collect::<Vec<_>>(),
            vec![(85, SectionType::Values), (145, SectionType::Raw)]
        );
    }
}
