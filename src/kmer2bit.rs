pub mod utils;

#[derive(Eq, Hash, PartialEq, Debug)]
pub struct Kmer2bit<const K: usize> {
    vec: [u8; K],
} 

impl<const K: usize> Kmer2bit<K> {

    /* 
    fn new() -> Self {
        Kmer2bit {
            vec: [0; K],
        }
    }
    */

    pub fn from_string(kmer: &str) -> Option<Self> {
    
        let mut twobit_array: [u8; K] = [0; K];

        for i in 0..K {
            let tkmer = &kmer[(4 * i)..(4 * i + 4)];
            let dna_quartet = utils::dna_quartet2u8(&tkmer);
            match dna_quartet {
                Some(n) => twobit_array[i] = n,
                _ => return None, 
            }
        
            // twobit_array[i] = utils::dna_quartet2u8(&tkmer)
        }
        Some(Kmer2bit{vec: twobit_array})
    }

    pub fn get_kmer_string(&self) -> Option<String> {
 
        let mut dna_vec: Vec<String> = Vec::new();
        for i in 0..K {
            dna_vec.push(utils::u82dna_quartet(self.vec[i]).expect("Input num out of the range!").to_string());
        }
        Some(dna_vec.join(""))
    }
}

