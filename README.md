

<img src="https://github.com/xjtu-omics/Mako/blob/master/supports/mako_logo.png" alt="mako_logo" width="30%" height="30%" align=center/>


Mako is a bottom-up guided model-free CSV detection tool. It first builds a mutational signal graph and utilizes pattern growth to detect maximal subgraphs as CSVs.

<img src="https://github.com/xjtu-omics/Mako/blob/master/supports/Mako_workflow.png" alt="mako_workflow" width="80%" height="80%" align=center/>

Please check the [wiki](https://github.com/xjtu-omics/Mako/wiki) page for more details.

## License

SVision is free for non-commercial use by academic, government, and non-profit/not-for-profit institutions. A commercial version of the software is available and licensed through Xi’an Jiaotong University.
For more information, please contact with Jiadong Lin (jiadong324@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn).

## Install and run

Mako requires Java JDK (>=1.8), we provide a prebuilt JAR package **Mako.jar** for directly usage. 
Please check [release](https://github.com/xjtu-omics/Mako/releases).

### Dependency

- htsjdk (https://github.com/samtools/htsjdk): A Java API for processing high-throughput sequencing (HTS) data.
- Python (V>=3.6): This is required for creating Mako configuration file. 
  - Required package: pysam, pandas, numpy

### Usage

**NOTE:** BAM file should under your working directory.
```
# Configuration
python ParseMako.py config -b sample.bam -n 30000 -w ./working_dir/ -s sampleName -f /path/to/ref.fa.fai
```

```
# Detection
java -jar Mako.jar -R /path/to/ref.fa -F /path/to/sampleName.mako.cfg
```

```
# Convert to VCF format (optional)
python ParseMako.py tovcf -m sampleName_mako_calls.txt -o sampleName_mako.vcf -r /path/to/ref.fa -s sampleName
```

### Run demo data

```
# Create configuration file
python ParseMako.py config -b NA19240.30X.chr20.1000K-2000K.bam -n 30000 -w ./working_dir/ -s NA19240 -f /demo.fa.fai

# Run Mako
java -jar /path/to/Mako.jar -R /path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa -F /path/to/NA19240.mako.cfg
```

## Known issues

1. Please make sure the reference used for running Mako is identical to the alignment one.
2. ...


## Citation
The manuscript is accepted by Genomics Proteomics Bioinformatics, which will be online soon.

## Contact

If you have any questions, please feel free to contact with Jiadong Lin (jiadong66@stu.xjtu.edu.cn) or Kai Ye (kaiye@xjtu.edu.cn)
