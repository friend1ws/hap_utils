use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::collections::{HashMap, HashSet};
use std::str;
use bio::io::{fasta, fastq};
use bio::alphabets::dna::revcomp;
use clap::{Parser, Subcommand};
use csv;

mod kmer2bit;

const BIT_SIZE: usize = 6;


#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
#[clap(propagate_version = true)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Extract haplotype specific kmer
    Parse {
        #[clap(value_parser)]
        hap1_file_path: String,
        hap2_file_path: String,
        out_file_path: String,        
    },
    Binning {
        #[clap(value_parser)]
        fasta_file_path: String,
        kmer_file_path: String,
        out_file_path: String,
    },
}

fn parse_main(hap1_file_path: &String, hap2_file_path: &String, out_file_path: &String) -> Result<(), Box<dyn Error>> {

    let mut kmer2count: HashMap<kmer2bit::Kmer2bit::<BIT_SIZE>, [u32; 2]> = HashMap::new();

    let reader = fasta::Reader::from_file(hap1_file_path)?;
    for result in reader.records() {
        let record = result?;
    
        let tseq_str = String::from_utf8(record.seq().to_vec())?;
        let seq_len = tseq_str.len();

        eprintln!("Processing {} {}", record.id(), seq_len);
        for i in 0..(seq_len - BIT_SIZE * 4) {
            let tkmer = &tseq_str[i..(i + BIT_SIZE * 4)];
            let trkmer: String = String::from_utf8(revcomp(tkmer.as_bytes())).unwrap();
            let trkmer: &str = &trkmer[..];

            match kmer2bit::Kmer2bit::<BIT_SIZE>::from_string(tkmer) {
                Some(tkmer2bit) => {
                    let count = (kmer2count).entry(tkmer2bit).or_insert([0, 0]);
                    (*count)[0] += 1;
                },
                None => continue,
            }

            match kmer2bit::Kmer2bit::<BIT_SIZE>::from_string(trkmer) {
                Some(tkmer2bit) => {
                    let count = (kmer2count).entry(tkmer2bit).or_insert([0, 0]);
                    (*count)[0] += 1;
                },
                None => continue,
            }

        }

    }

    let reader = fasta::Reader::from_file(hap2_file_path)?;
    for result in reader.records() {
        let record = result?;

        let tseq_str = String::from_utf8(record.seq().to_vec())?;
        let seq_len = tseq_str.len();

        eprintln!("Processing {} {}", record.id(), seq_len);
        for i in 0..(seq_len - BIT_SIZE * 4) {
            let tkmer = &tseq_str[i..(i + BIT_SIZE * 4)];

            let trkmer: String = String::from_utf8(revcomp(tkmer.as_bytes())).unwrap();
            let trkmer: &str = &trkmer[..];

            match kmer2bit::Kmer2bit::<BIT_SIZE>::from_string(tkmer) {
                Some(tkmer2bit) => {
                    let count = (kmer2count).entry(tkmer2bit).or_insert([0, 0]);
                    (*count)[1] += 1;
                },
                None => continue,
            }
            match kmer2bit::Kmer2bit::<BIT_SIZE>::from_string(trkmer) {
                Some(tkmer2bit) => {
                    let count = (kmer2count).entry(tkmer2bit).or_insert([0, 0]);
                    (*count)[1] += 1;
                },
                None => continue,
            }

        }

    }

    let write_file = File::create(out_file_path).unwrap();
    let mut writer = BufWriter::new(&write_file);


    for (kmer, count) in &kmer2count {
        let kmer_string = kmer.get_kmer_string().unwrap();
        if ( (*count)[0] == 0 && (*count)[1] > 0 ) || ((*count)[0] > 0 && (*count)[1] == 0 ) {
            writeln!(&mut writer, "{}\t{}\t{}", kmer_string, (*count)[0], (*count)[1])?;
        }
    }

    Ok(())

}


fn binning_main(fasta_file_path: &String, kmer_file_path: &String, out_file_path: &String) -> Result<(), Box<dyn Error>> {

    let mut kmer_hap1: HashSet<String> = HashSet::new();
    let mut kmer_hap2: HashSet<String> = HashSet::new();

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(File::open(kmer_file_path)?);

    for result in rdr.records() {
        let record = result?;
        let kmer: String = record[0].to_string();
        let count1: u32 = record[1].parse()?;
        let count2: u32 = record[2].parse()?;

        if count1 > 0 {
            kmer_hap1.insert(kmer);
        } else if count2 > 0 {
            kmer_hap2.insert(kmer);
        }
    }
   

    let reader = fastq::Reader::from_file(fasta_file_path)?;
    let write_file = File::create(out_file_path).unwrap();
    let mut writer = BufWriter::new(&write_file);

    for result in reader.records() {
        let record = result?;

        let tseq_str = String::from_utf8(record.seq().to_vec())?;
        let seq_len = tseq_str.len();

        let mut hap_kmer_count = vec![0, 0];

        for i in 0..(seq_len - BIT_SIZE * 4) {
            let tkmer = &tseq_str[i..(i + BIT_SIZE * 4)];
            
            let tkmer_string: String = tkmer.to_string();
            if kmer_hap1.contains(&tkmer_string) {
                hap_kmer_count[0] = hap_kmer_count[0] + 1
            }
            if kmer_hap2.contains(&tkmer_string) {
                hap_kmer_count[1] = hap_kmer_count[1] + 1
            }

        }

        writeln!(&mut writer, "{}\t{}\t{}", record.id(), hap_kmer_count[0], hap_kmer_count[1])?;

    }

    Ok(())

}

fn main() -> Result<(), Box<dyn Error>> {

    let cli = Cli::parse();

    let result;
    match &cli.command {
        Commands::Parse { hap1_file_path, hap2_file_path, out_file_path } => {
            result = parse_main(hap1_file_path, hap2_file_path, out_file_path);
        },
        Commands::Binning { fasta_file_path, kmer_file_path, out_file_path } => {
            result = binning_main(fasta_file_path, kmer_file_path, out_file_path);
        },
    }
    
    result
}


