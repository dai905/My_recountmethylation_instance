[![DOI](https://zenodo.org/badge/337589428.svg)](https://zenodo.org/badge/latestdoi/337589428)

# recountmethylation_instance

Authors: Sean Maden, Reid Thompson, Kasper Hansen, Abhi Nellore

Set up and maintain an instance, or synchronization, of public DNA methylation (DNAm) arrays from the 
[Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/). The `recountmethylation_instance` 
resource provides a `snakemake` workflow including all steps to create HDF5 `.h5` and HDF5-SummarizedExperiment 
`h5se` database files of compiled and harmonized DNAm array assays and sample metadata. Expedited setup of a 
synchronization environment is supported using `sh` or `conda` (Tutorial, below).

# Who is this resource for?

The `recountmethylation_instance` resource can handle data compilations comprising upwards of tens of 
thousands of samples and hundreds of studies. If you are instead obtaining data from just one or a few 
samples or studies from GEO, there are several more lightweight approaches available to you, including 
manual download of files from the GEO website, and programmatic options including functions provided in the [`GEOquery`](https://www.bioconductor.org/packages/release/bioc/html/GEOquery.html) library, the `getGenomicRatioSetFromGEO()` 
function from the [`minfi`](http://bioconductor.org/packages/release/bioc/html/minfi.html) library, and 
the `gds_idatquery()` and `gds_idat2rg()` functions from the [`recountmethylation`](https://www.bioconductor.org/packages/devel/bioc/html/recountmethylation.html) library.

# Dependencies

Dependencies to run an instance are shown in the below tables showing dependency name (column 1), type (column 2), 
confirmed version (column 3), and requirement (column 4). Instances were confirmed to work with the indicated 
versions (column 2), and in many cases can work with more recent versions. Several dependencies are optional to 
run specific features of an instance, as indicated in the "Required?" column.

Main software and programming language dependencies:

|Dependency        |         Type         |  Version  | Required? |
|------------------|----------------------|-----------|-----------|
|  [Python 3](https://www.python.org/downloads/)      | programming language |  >=3.7.3  |    yes    |
|     [R](https://cran.r-project.org/)        | programming language |  >=3.6.0  |    yes    |
| [RabbitMQ](https://www.rabbitmq.com/)     |        broker        |   3.8.11  |    yes    |
|  [MongoDB](https://www.mongodb.com/2)     |      db syntax       |   4.4.3   |    yes    |
|  [SQLite](https://www.sqlite.org/download.html)      |      db syntax       |   3.30.1  |    yes    |
| [Entrez Utilities](https://www.ncbi.nlm.nih.gov/books/NBK25500/) | NCBI/GEO API utils   |   10.9    |    yes    |
|  [Python 2](https://www.python.org/downloads/)      | programming language |   2.7.5   |    no     |

The following Python 3 libraries are required or recommended. Note packages 
`dash` and `plotly` are only required to run the optional server dashboard utility:

|  Library  |  Type  | Version | Required? |
|-----------|----|---------|-----------|
| [snakemake](https://pypi.org/project/snakemake/) |  Python 3 library |  6.1.2  |  yes  |
|   [pandas](https://pypi.org/project/pandas/)  |  Python 3 library |  0.25.1 |  yes  |
|   [celery](https://pypi.org/project/celery/)  |  Python 3 library |  4.2.1  |  yes  |
|  [pymongo](https://pypi.org/project/pymongo/)  |  Python 3 library |  3.7.2  |  yes  |
|    [dash](https://pypi.org/project/dash/)   |  Python 3 library |  1.20.0 |  no  |
|   [plotly](https://pypi.org/project/plotly/)  |  Python 3 library |  4.14.3 |  no  |

The following R libraries are required or recommended. Note the 
`ewastools` library is only required for optional genotype-based shared 
identity analyses, which are made as part of the calculated DNAm-based metadata.

|   Library   |  Type  | Version |  Required?  |
|-------------|--------|---------|-------------|
|    [minfi](https://bioconductor.org/packages/release/bioc/html/minfi.html) | R library |  1.32.0 |    yes   |
|    [rhdf5](https://www.bioconductor.org/packages/release/bioc/html/rhdf5.html)    | R library |  2.30.1 |   yes |
| [DelayedArray](https://www.bioconductor.org/packages/release/bioc/html/DelayedArray.html)| R library |  0.12.3 |   yes |
|  [HDF5Array](https://www.bioconductor.org/packages/release/bioc/html/HDF5Array.html)  | R library |  1.14.4 |  yes |
|  [ewastools]()  |  R library |   1.7   |  no  |

Python 2 is only required if you intend to run the `MetaSRA-pipeline`. For this, 
the following libraries are recommended:

|    Library   |  Type   | Version   | Required? |
|--------------|---------|-----------|-----------|
|    [numpy](https://pypi.org/project/numpy/)     |  Python 2 library  | 1.15.4  |  no  |
|    [scipy](https://pypi.org/project/scipy/)     |  Python 2 library  | 1.1.0  |  no  |
| [scikit-learn](https://pypi.org/project/scikit-learn/) |  Python 2 library  | 0.20.1  |  no  |
|  [setuptools](https://pypi.org/project/setuptools/)  |  Python 2 library  |  0.9.8  |  no  |
|  [marisa-trie](https://pypi.org/project/marisa-trie/) |  Python 2 library  |  0.7.5  |  no  |
|     [dill](https://pypi.org/project/dill/)    |  Python 2 library  | 0.2.8.2 |  no  |
|     [nltk](https://pypi.org/project/nltk/)     |  Python 2 library  |   3.4   |  no  |

# Setup

First clone the latest version of the `recountmethylation_instance` repo from GitHub, then 
navigate into the directory.

```
git clone https://github.com/metamaden/recountmethylation_instance
cd recountmethylation_instance
```

Next, choose a setup option below. Several options have been provided to support various 
OS environments.

## 1a. Setup with sh script

For setup with `sh`, run the provided script `setup_instance.sh`:

```
sh setup_instance.sh
```

## 1b. Setup with a conda virtual environment

For setup using virtual environments with `conda`, you may run the script `anaconda_setup.sh`:

```
sh anaconda_setup.sh
```

Alternatively, create the environment from the provided `.yml` file:

```
conda env create -f ./inst/yml/environment_rmi_py3.yml
```

This will set up the main environment, called `py3,` to run the core data processing rules. 

For running MetaSRA-pipeline, you may either set up an additional Python 2 environment or 
run the corresponding `snakemake` rules manually if a Python 2 alias is callable from the 
`py3` environment.

# Tutorial

This tutorial shows how to set up and initiate synchronization of public DNA methylation array data. 
Note that most steps involve calling `snakemake`rules defined in the provided `Snakefile` script. Logs, 
including stdout, for each rule called are stored by default in the "snakemakelogs" subdirectory.

## 1. Set up the instance

For this example, DNAm array data run using the Illumina Infinium HumanMethylation 450K array 
platform will be considered. First clone the resource from GitHub, then rename the cloned repo.

```
git clone https://github.com/metamaden/recountmethylation_instance
mv recountmethylation_instance recountmethylation_instance_hm450k
cd recountmethylation_instance_hm450k
```

Next, we can perform setup with `setup_instance.sh`:

```
sh setup_instance.sh
```

After successfully completing the above setup steps, the top level of the new instance should look something like this:

![Top level of an instance after setup.](img/example_rmi_toplvl.PNG)

## 2. Configure the instance
 
Next, we need to configure the instance, including specifying the array platform to target, and 
specifying sample IDs to exclude.

### 2a. Specify target platform

The HM450K platform is presently targeted by default, but this may change. To explicitly set
the platform to target, run the following:

```
python3 ./recountmethylation_server/src/set_acc.py
```

This produces the dialogue:

```
(Y/N) Current target platform accession ID is GPL13534.
Do you want to change the target platform?
```

Entering `Y` returns:

```
(1/2/3) Enter the number of the platform to target:
 1. HM450K (GPL13534)
 2. EPIC/HM850K (GPL21145)
 3. HM27K (GPL8490)
```

Type `1` to specify the HM450K platform.

### 2b. Run a new Entrez Utilities query

Data files are recognized by queries to the GEO DataSets API using Entrez Utilities software. Running 
a fresh query will identify all valid data files for the targeted platform. To do this, enter:

```
snakemake --cores 1 new_eqd
```

### 2c. Exclude freeze sample IDs

The sample/GSM IDs for the most recently available data freezes are included at `./inst/freeze_gsmv/`.
Excluding these GSM IDs for this instance will allow us to synchronize the subset of samples available 
since the latest data freeze (currently: November 7, 2020).

```
snakemake --cores 1 exclude_gsm
```

For this example instance, the rule generated the following output:

```
Starting with latest detected filter file: recount-methylation-files/equery/gsequery_filt.1618960590
Applying filter...
After filter, retained 72 studies.
Writing new filter file:  recount-methylation-files/equery/gsequery_filt.1619031416
```

Functions downstream will now recognize and use the newer version of the file `gsequery_filt.*` according 
to the newer applied timestamp `1619031416`.

### 2d. Initialize the instance metadata

Use the following rule to initialize the metadata, including timestamp and version, for this instance:

```
snakemake --cores 1 new_instance_md
```

This creates a new subdirectory containing the instance metadata. This is also where the newly
generated metadata files will be stored, including metadata mapped from GSM JSON files and 
DNAm model-based predicted metadata.

## 3. Running the server and downloading files

This section shows how to download datasets from the GEO servers to an instance. The files of interest 
include the SOFT files stored under a unique GSE accession for each study, and a pair of IDAT files 
for each unique GSM accession corresponding to a sample. These files are compressed, so you can expect
new downloads to all have a `.gz` extension.

For example, after successfully downloading a series of study SOFT files, the directory at 
`./recount-methylation-files/gse_soft/` might look like this:

![Compressed SOFT files](img/example_soft-compressed.PNG)

And the directory at `./recount-methylation-files/idats/` might look something like this:

![Compressed IDAT files](img/example_idats-compressed.PNG)

`recountmethylation_instance` uses the Entrez Direct queries in the previous section to make URLs
to download study and sample datasets from the GEO servers. Running the server is a safe way to 
ensure the correct files are downloaded and that downloaded files are tracked. If the server will
not run continuously, or repeated attempts to complete server runs fail, you can still download the
files using the provided `.R` scripts.

### 3a. Running the server

Before we can start downloading public data, we need to run the MongoDB service with sudo access. 
This can be done with either:

```
service mongod start
```

or

```
sudo service mongod start
```

Once MongoDB is running, we can initialize the server with:

```
snakemake --cores 1 run_server
```

If this rule hangs, you may alternatively call the script manually:

```
python3 ./recountmethylation_server/src/server.py
```

The `server.py` script process will systematically target and download study SOFT files and sample IDAT files, 
according to the contents of the filtered EDirect query files.

Note, you may need to restart the server process periodically if your connection is interrupted, the MongoDB service 
stops, etc. To avoid repeated hanging on corrupt or malformed files, target study/GSE ids are shuffled for
each `server.py` run.

As the server process runs, you may monitor its progress. Check the total study SOFT and sample IDAT files downloaded with:

```
ls -l recount-methylation-files/gse_soft/ | wc -l
ls -l recount-methylation-files/idats/ | wc -l
```

For added convenience, a server dashboard utility has been provided. This displays the instance files over time, and allows
you to track the addition of new files over time. Run the dashboard with:

```
snakemake --cores 1 server_dash
```

Entering the displayed IP address into a browser tab will display the server dashboard.

### 3b. Downloading SOFT and IDAT files with provided `.R` scripts

An alternative method of obtaining study and sample data from the GEO servers is to use the provided 
`.R` scripts. These scripts compare downloaded files to the full set of study/GSE and sample/GSM IDs, 
and they will try downloading data for remaining IDs. To download remaining IDATs, use:

```
Rscript ./inst/scripts/dl_idats.R
```

To download SOFT files, run:

```
Rscript ./inst/scripts/dl_soft.R
```

## 4. Reformatting and compiling data files

Once data files have been downloaded, they can be prepared for compilation. 

Over the course of the following steps, you may find unexpected file types have been downloaded, such as files with 
`.txt.gz` extension that were downloaded to the `./recount-methylation-files/idats` folder (see below). These types of 
anomalies are to be expected when performing a comprehensive acquisition and compilation of public data. 

![Expanded IDATs and .txt files.](img/example_idats-expand.PNG)

The instance will be able to handle many of these types of anomalies, however it is still possible new anomalous filenames and
extensions can be encountered whenever new public data is compiled. It is entirely up to you, the user, to decide how to handle 
files in these situations. For instance, you may choose to manually delete these files, or alternatively you could move them to 
a new location (e.g. move `.txt` files to a new directory at  `./recount-methylation-files/txt/`) and then target them in a future 
compilation.

### 4a. Sample IDATs

The sample IDATs are paired files containing red and green signals from the array runs. Thus, two 
valid IDATs are expected per sample/GSM id, where their filenames end as `.*_Red.idat.*` and `.*_Grn.idat.*`.

Since we download the IDATs as compressed `.gz` files, we need to expand them.

```
snakemake --cores 1 unzip_idats
```

Because paired IDATs may have distinct filename timestamps, it is necessary to create hardlinks to ensure 
timestamps match for sample IDAT pairs.

```
snakemake --cores 1 make_idat_hlinks
```

Once the IDATs have been prepared, they can be compiled into various database file types. To begin,
read the raw/unnormalized red and green signals are into flat compilation tables with `get_rg_compilations`. 

```
snakemake --cores 1 get_rg_compilations
```

Since downloaded IDATs may occasionally be malformed, mislabeled, or otherwise unreadable, a series of checks
are performed automatically, including an evaluation of file sizes and similar sizes for paired red and 
green signal files. In practice, this means an initial attempt is made to read in samples as batches, or chunks. 
If the attempt to read a given chunk fails, individual sample IDATs are read in successively and combined, and 
invalid files are excluded.

Next, the flat red and green signal compilation tables are used to make HDF5 `h5` and HDF5-SummarizedExperiment 
`h5se` objects.

```
snakemake --cores 1 get_h5db_rg
snakemake --cores 1 get_h5se_rg
```

Once the red and green `RGChannelSet` tables and database files are generated, you may proceed to generate 
additional compilation files. The full list of compilation rules, in their recommended order, is as follows:

1. `get_rg_compilations`: Make flat tables containing the Red and Green signal intensities. 
2. `get_h5db_rg`: Make an HDF5 database `.h5` file containing the Red and Green signal intensities. 
3. `get_h5se_rg`: Make an HDF5-SummarizedExperiment `h5se` file containing the Red and Green signal intensities
4. `get_h5db_gm`: Make an HDF5 database `.h5` file containing the Methylated and Unmethylated signals.
5. `get_h5se_gm`: Make an HDF5-SummarizedExperiment `h5se` file containing the Methylated and Unmethylated signals.
6. `get_h5db_gr`: Make an HDF5 database `.h5` file containing the Beta-values (DNAm fractions).
7. `get_h5se_gr`: Make an HDF5-SummarizedExperiment `h5se` file containing the Beta-values (DNAm fractions).

You could run the above from an active R session by running `R` and then entering the following:

```
library(recountmethylation.pipeline)
get_rg_dtables()
get_h5db_rg()
get_h5se_rg()
get_h5db_gm()
get_h5se_gm()
get_h5db_gr()
get_h5se_gr()
```

You could also run the following one-liner from your Terminal to call the steps in their recommended order:

```
snakemake --cores 1 run_dnam_pipeline
```

### 4b. Study SOFT files

Sample metadata is contained in the SOFT files. After expanding the `.gz` compressed SOFT files, we need to 
extract the sample-specific metadata into `.json` files before mapping with either MetaSRA-pipeline or the 
included mapping scripts. After extraction, the `.json` files are further filtered to remove study-specific 
metadata.

```
snakemake --cores 1 process_soft
snakemake --cores 1 apply_jsonfilt
```

As an example, SOFT-extracted metadata for [`GSM4671807`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4671807) 
would be stored at `./gsm_soft/` and look something like:

```
!Sample_title = K562 - control 1

!Sample_geo_accession = GSM4671807

!Sample_status = Public on Nov 23 2022

...
...

!Sample_characteristics_ch1 = cell line: K562

!Sample_characteristics_ch1 = disease: chronic myelogenous leukemia

!Sample_characteristics_ch1 = treatment: control
```

By contrast, the JSON-formatted version of `GSM4671807`'s metadata would be stored at `./gsm_json/` and appear as:

```
[
  {
    "!Sample_title": "K562 - control 1",
    "!Sample_geo_accession": "GSM4671807",
    "!Sample_status": "Public on Nov 23 2022",
    "!Sample_submission_date": "Jul 15 2020",
    "!Sample_last_update_date": "Nov 23 2022",
    "!Sample_type": "genomic",
    "!Sample_channel_count": "1",
    "!Sample_source_name_ch1": "control replicate 1",
    "!Sample_organism_ch1": "Homo sapiens",
    "!Sample_taxid_ch1": "9606",
    "!Sample_characteristics_ch1": "cell line: K562",
    "!Sample_characteristics_ch1.1": "disease: chronic myelogenous leukemia",
    "!Sample_characteristics_ch1.2": "treatment: control",
    "!Sample_treatment_protocol_ch1": "GSE154471\n",
    "cell_line": "K562",
    "treatment": "control"
  }
]
```

And finally the filtered JSON file for `GSM4671807`, stored at `./gsm_jsonfilt/`, looks like:

```{r}
[
{
  "!Sample_characteristics_ch1": "cell line: K562",
  "!Sample_characteristics_ch1.1": "disease: chronic myelogenous leukemia",
  "!Sample_characteristics_ch1.2": "treatment: control",
  "!Sample_source_name_ch1": "control replicate 1",
  "!Sample_title": "K562 - control 1"
}
]
```

This JSON formatting is important because it allows the file to be read by the `MetaSRA-pipeline`. We also
apply a filter on metadata fields such as `!Sample_treatment_protocol_ch1` which are likely to be repeated 
across samples. This can help ensure the pipeline maps characteristics which differentiate samples within
the same experiment, rather than mapping the same terms to every sample in a given experiment.

Once the `.json` files with sample-specific metadata have been prepared, you have the option of running
any of the following available metadata processing rules:

* `do_mdmap`: Map and harmonize metadata using the provided scripts. These scripts use regular expressions to 
            automatically detect and categorize tags in `.json` files, and then to uniformly format and 
            annotate metadata terms under several columns, including "disease" (e.g. disease condition or 
            experiment group), "tissue" (e.g. tissue of origin), "age" (chronological age), and "sex" 
            (provided sex information).
* `run_msrap`: Run the MetaSRA-pipeline. This produces sample type predictions, as well as ENCODE ontology terms 
             from several major ontology dictionaries.
* `do_dnam_md`: Get DNAm-derived metadata (model-based predictions for age, sex, and blood cell types) and quality 
              metrics (BeadArray controls, methylated signal, unmethylated signal, and predictated replicates). 
              This should be run after all 3 types of h5se compilations are complete.
              
Once one or all of these rules have been successfully run, compile and append the harmonized metadata to the available DNAm data compilations:

```
snakemake --cores 1 make_md_final
snakemake --cores 1 append_md
```

You could also accomplish the above by running `R` and then entering the following:

```{r}
library(recountmethylation.pipeline)
get_mdmap()
get_mddnam()
md_agg()
append_md_snakemake()
```

# Additional resources

If you find `recountmethylation_instance` useful, the following resources may also be helpful.

## Packages

* `recountmethylation`. R/Bioconductor package providing utilities to access and analyze public DNAm array data from GEO. [10.18129/B9.bioc.recountmethylation](https://doi.org/doi:10.18129/B9.bioc.recountmethylation)

* `minfi`. R/Bioconductor package for analysis of DNA methylation arrays. [10.18129/B9.bioc.minfi](10.18129/B9.bioc.minfi )

* `HDF5Array`. R/Bioconductor package for handling HDF5 and HDF5SummarizedExperiment objects with DelayedArray-powered backends. [10.18129/B9.bioc.HDF5Array](10.18129/B9.bioc.HDF5Array) 

## Papers

* Sean K Maden, Reid F Thompson, Kasper D Hansen, Abhinav Nellore, *Human methylome variation across Infinium 450K data on the Gene Expression Omnibus*, NAR Genomics and Bioinformatics, Volume 3, Issue 2, June 2021, https://doi.org/10.1093/nargab/lqab025 

* Timothy J Triche Jr, Daniel J Weisenberger, David Van Den Berg, Peter W Laird, Kimberly D Siegmund, *Low-level processing of Illumina Infinium DNA Methylation BeadArrays*, Nucleic Acids Research, Volume 41, Issue 7, April 2013, https://doi.org/10.1093/nar/gkt090
