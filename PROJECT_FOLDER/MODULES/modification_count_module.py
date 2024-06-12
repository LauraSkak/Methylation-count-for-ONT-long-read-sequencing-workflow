
from gwf import AnonymousTarget

def run_modkit(alignment_file, modkit_outfile, reference):
    """
    Takes an unphased .bam file and creates a methylation count table.
    """
    inputs = [alignment_file]
    outputs = [f'{modkit_outfile}.gz', f'{modkit_outfile}.gz.tbi']
    options = {"walltime":"10:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {outfile} --ref {ref} --threads 16"

    modkit pileup \
        {infile} {outfile} \
        --ref {ref} \
        --only-tabs \
        --threads 16
    
    conda activate samtools
    
    echo "bgzip {outfile}"
    
    bgzip {outfile} 
    
    echo "tabix -p bed -f {outfile}.gz "
    
    tabix -p bed -f {outfile}.gz 
    
    '''.format(infile = alignment_file, outfile = modkit_outfile, ref = reference)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def run_modkit_merged_strand(alignment_file, modkit_outfile, reference):
    """
    Takes an unphased .bam file and creates a methylation count table, where strands are merged. The table only output modifications placed on a CpG site.
    """
    inputs = [alignment_file]
    outputs = [f'{modkit_outfile}.gz', f'{modkit_outfile}.gz.tbi']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {outfile} --ref {ref} --threads 16 --combine-strands --only-tabs --cpg"

    modkit pileup \
    {infile} {outfile} \
    --ref {ref} \
    --threads 16 \
    --combine-strands \
    --only-tabs \
    --cpg

    conda activate samtools
    
    echo "bgzip {outfile}"
    
    bgzip {outfile} 
    
    echo "tabix -p bed -f {outfile}.gz "
    
    tabix -p bed -f {outfile}.gz 
    
    '''.format(infile = alignment_file, outfile = modkit_outfile, ref = reference)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def run_modkit_merged_strand_and_mods(alignment_file, modkit_outfile, reference):
    """
    Takes an unphased .bam file and creates a methylation count table, where strand and modification types are merged. The table only output modifications placed on a CpG site.
    """
    inputs = [alignment_file]
    outputs = [f'{modkit_outfile}.gz', f'{modkit_outfile}.gz.tbi']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {outfile} --ref {ref} --threads 16 --combine-mods --combine-strands --only-tabs --cpg"

    modkit pileup \
    {infile} {outfile} \
    --ref {ref} \
    --threads 16 \
    --combine-mods \
    --combine-strands \
    --only-tabs \
    --cpg

    conda activate samtools
    
    echo "bgzip {outfile}"
    
    bgzip {outfile} 
    
    echo "tabix -p bed -f {outfile}.gz "
    
    tabix -p bed -f {outfile}.gz 
    
    '''.format(infile = alignment_file, outfile = modkit_outfile, ref = reference)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def run_modkit_on_regions(alignment_file, modkit_outfile, region_file, reference):
    """
    Takes an unphased .bam file and creates a methylation count table, where strand and modification types are merged. The table only output modifications placed on a CpG site.
    """
    inputs = [alignment_file]
    outputs = [f'{modkit_outfile}.gz', f'{modkit_outfile}.gz.tbi']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {outfile} --ref {ref} --include-bed {cpgs} --threads 16 --combine-mods --combine-strands --only-tabs --cpg"

    modkit pileup \
    {infile} {outfile} \
    --ref {ref} \
    --include-bed {cpgs} \
    --threads 16 \
    --only-tabs \
    --cpg

    conda activate samtools
    
    echo "bgzip {outfile}"
    
    bgzip {outfile} 
    
    echo "tabix -p bed -f {outfile}.gz "
    
    tabix -p bed -f {outfile}.gz 
    
    '''.format(infile = alignment_file, outfile = modkit_outfile, ref = reference, cpgs = region_file)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def run_modkit_on_regions_merged_strand(alignment_file, modkit_outfile, region_file, reference):
    """
    Takes an unphased .bam file and creates a methylation count table, where strand and modification types are merged. The table only output modifications placed on a CpG site.
    """
    inputs = [alignment_file]
    outputs = [f'{modkit_outfile}.gz', f'{modkit_outfile}.gz.tbi']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {outfile} --ref {ref} --include-bed {cpgs} --threads 16 --combine-mods --combine-strands --only-tabs --cpg"

    modkit pileup \
    {infile} {outfile} \
    --ref {ref} \
    --include-bed {cpgs} \
    --threads 16 \
    --combine-strands \
    --only-tabs \
    --cpg

    conda activate samtools
    
    echo "bgzip {outfile}"
    
    bgzip {outfile} 
    
    echo "tabix -p bed -f {outfile}.gz "
    
    tabix -p bed -f {outfile}.gz 
    
    '''.format(infile = alignment_file, outfile = modkit_outfile, ref = reference, cpgs = region_file)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def run_modkit_on_regions_merged_strand_and_mods(alignment_file, modkit_outfile, region_file, reference):
    """
    Takes an unphased .bam file and creates a methylation count table, where strand and modification types are merged. The table only output modifications placed on a CpG site.
    """
    inputs = [alignment_file]
    outputs = [f'{modkit_outfile}.gz', f'{modkit_outfile}.gz.tbi']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 16}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {outfile} --ref {ref} --include-bed {cpgs} --threads 16 --combine-mods --combine-strands --only-tabs --cpg"

    modkit pileup \
    {infile} {outfile} \
    --ref {ref} \
    --include-bed {cpgs} \
    --threads 16 \
    --combine-mods \
    --combine-strands \
    --only-tabs \
    --cpg

    conda activate samtools
    
    echo "bgzip {outfile}"
    
    bgzip {outfile} 
    
    echo "tabix -p bed -f {outfile}.gz "
    
    tabix -p bed -f {outfile}.gz 
    
    '''.format(infile = alignment_file, outfile = modkit_outfile, ref = reference, cpgs = region_file)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def phasing_modkit(alignment_file, reference, prefix, outfile_dir):
    """
    Takes a phased .bam file and creates a methylation count table for each haplotype.
    """
    inputs = [alignment_file, f'{alignment_file}.bai']
    outputs = [f'{outfile_dir}/{prefix}_1.bed', f'{outfile_dir}/{prefix}_2.bed', f'{outfile_dir}/{prefix}_ungrouped.bed', f'{outfile_dir}/{prefix}.DONE']
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb", "cores": 32}
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate ONT

    echo "Job ID: $SLURM_JOB_ID"
    echo "modkit pileup {infile} {prefix}.bed --ref {ref} --partition-tag HP --prefix {prefix} --threads 32"

    modkit pileup \
        {infile} \
        {out_dir} \
        --ref {ref} \
        --partition-tag HP \
        --prefix {prefix} \
        --combine-strands \
        --cpg \
        --only-tabs \
        --threads 32
        
    touch {out_dir}/{prefix}.DONE
    
    '''.format(ref = reference, infile = alignment_file, prefix=prefix, out_dir = outfile_dir)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def merge_modkit(infile_template, outfile, regions_bed_file, samples_meta_file, keep_sites):
    """
    FIXME
    """
    inputs = [] #FIXME: Missing infline
    outputs = [outfile]
    options = {"walltime":"12:00:00","account":"sexChromosomes", "memory":"50gb"}
    
    if keep_sites:
        keep_sites = "--keep_sites "
    else:
        keep_sites = ""
    
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base

    echo "Job ID: $SLURM_JOB_ID"
    echo "python aggregate_site_data_per_region.py \
        --infile {infile} \
        --outfile {outfile} \
        --regions {bed_file} \
        --sample_meta {meta_file} \
        --keep_sample_info \
        --keep_mod \
        {keep_sites}"

    cd SCRIPTS
    
    python aggregate_site_data_per_region.py \
        --infile {infile} \
        --outfile {outfile} \
        --regions {bed_file} \
        --sample_meta {meta_file} \
        --keep_sample_info \
        --keep_mod \
        {keep_sites}
    
    '''.format(infile = infile_template, outfile = outfile, bed_file = regions_bed_file, meta_file = samples_meta_file, keep_sites = keep_sites)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def merge_phased_modkit(infile_list, infile_template, outfile, regions_bed_file, samples_meta_file, keep_sites):
    """
    FIXME
    """
    inputs = infile_list
    outputs = [outfile]
    options = {"walltime":"24:00:00","account":"sexChromosomes", "memory":"50gb"}
    
    if keep_sites:
        keep_sites = "--keep_sites "
    else:
        keep_sites = ""
        
    spec = '''
    
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate base

    echo "Job ID: $SLURM_JOB_ID"
    echo "python aggregate_site_data_per_region.py \
        --infile {infile} \
        --outfile {outfile} \
        --regions {bed_file} \
        --sample_meta {meta_file} \
        --keep_sample_info \
        --phased \
        --keep_mod \
        {keep_sites}"
    
    cd SCRIPTS
    
    python aggregate_site_data_per_region.py \
        --infile {infile} \
        --outfile {outfile} \
        --regions {bed_file} \
        --sample_meta {meta_file} \
        --keep_sample_info \
        --phased \
        --keep_mod \
        {keep_sites}
    
    '''.format(infile = infile_template, outfile = outfile, bed_file = regions_bed_file, meta_file = samples_meta_file, keep_sites = keep_sites)
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)