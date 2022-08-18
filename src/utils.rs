//! Utils function for KFF

/* std use */

/* crate use */

/* project use */

/* mod declaration */
pub mod global_index;
pub mod kmer;
pub mod read;
pub mod write;

/* pub use */
pub use global_index::GlobalIndex;
pub use kmer::{Data, Kmer, Seq2Bit};
pub use read::KffRead;
pub use write::KffWrite;
#[inline]
pub(crate) fn ceil_to_8(n: u64) -> u64 {
    (n + 7) & !(7)
}

#[inline]
pub(crate) fn bits2store_k(k: u64) -> u64 {
    k * 2
}

#[inline]
pub(crate) fn bytes2store_k(k: u64) -> u64 {
    ceil_to_8(bits2store_k(k)) / 8
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::error;
    use rand::distributions::Distribution;

    #[test]
    fn ceil_to_8_() -> error::Result<()> {
        assert_eq!(ceil_to_8(1), 8);
        assert_eq!(ceil_to_8(7), 8);

        assert_eq!(ceil_to_8(8), 8);

        assert_eq!(ceil_to_8(9), 16);
        assert_eq!(ceil_to_8(15), 16);

        Ok(())
    }

    #[test]
    fn bits2store_k_() -> error::Result<()> {
        let range = rand::distributions::Uniform::from(0..100);
        let mut rng = rand::thread_rng();

        for _ in 0..100 {
            let value = range.sample(&mut rng);

            assert_eq!(bits2store_k(value), value * 2);
        }

        Ok(())
    }

    #[test]
    fn bytes2store_k_() -> error::Result<()> {
        assert_eq!(bytes2store_k(1), 1);
        assert_eq!(bytes2store_k(4), 1);
        assert_eq!(bytes2store_k(5), 2);
        assert_eq!(bytes2store_k(16), 4);
        assert_eq!(bytes2store_k(17), 5);

        Ok(())
    }
}
