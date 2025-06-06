process pystats {
    input:
        val mypath
    output:
        stdout
        //val mypath
        //path "pyoutputs.txt", emit: pyoutputs
        
    $/
    #!/usr/bin/env python3
    import subprocess
    
    items = "${mypath}".strip().split("/")
    #print(items[-1])
    filepath1 = "${mypath}"+"/alignment/"+items[-1]+".coverage.txt"
    #print(filepath1)
    with open(filepath1, 'r') as cov_report:
        header = cov_report.readline()
        header = header.rstrip()
        stats = cov_report.readline()
        stats = stats.rstrip()
        stats = stats.split()
        ref_name = stats[0]
        #print(ref_name)
        start = stats[1]
        end = stats[2]
        reads_mapped = stats[3]
        cov_bases = stats[4]
        cov = stats[5]
        depth = stats[6]
        baseq = stats[7]
        #print(reads_mapped)
        mapq = stats[8]
        
    #Get number of raw reads
    proc_1 = subprocess.run('zcat ' + "${mypath}/" + items[-1] + '_1_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_1 = proc_1.stdout.rstrip()
    reads_1 = int(wc_out_1) / 4
    proc_2 = subprocess.run('zcat ' + "${mypath}/" + items[-1] + '_2_humanclean.fastq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_2 = proc_2.stdout.rstrip()
    reads_2 = int(wc_out_2) / 4
    raw_reads = reads_1 + reads_2
    raw_reads = int(raw_reads)

    #Get number of clean reads
    proc_c1x = subprocess.run('zcat ' + "${mypath}/" + items[-1] + '_1.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_c1x = proc_c1x.stdout.rstrip()
    reads_c1x = int(wc_out_c1x) / 4
    proc_c2x = subprocess.run('zcat ' + "${mypath}/" + items[-1] + '_2.fq.gz | wc -l', shell=True, capture_output=True, text=True, check=True)
    wc_out_c2x = proc_c2x.stdout.rstrip()
    reads_c2x = int(wc_out_c2x) / 4
    clean_reads = reads_c1x + reads_c2x
    clean_reads = int(clean_reads)
    #print(clean_reads)
    
    #Get percentage of mapped reads/clean reads
    percent_map = "%0.4f"%((int(reads_mapped)/int(clean_reads))*100)
    #print(percent_map)
    
    #Gather QC metrics for consensus assembly
    filepath2 = "${mypath}"+"/assembly/"+items[-1]+".consensus.fa"
    with open(filepath2, 'r') as assem:
        header = assem.readline()
        header = header.rstrip()
        bases = assem.readline()
        bases = bases.rstrip()
        num_bases = len(bases)
        ns = bases.count('N')
        called = num_bases - ns
        pg = "%0.4f"%((called/int(end))*100)
        #print(called)
        #print(end)
        #print(pg)
    #Rename header in fasta to just sample name
    subprocess.run("sed -i \'s/^>.*/>"+items[-1]+"/\' "+filepath2, shell=True, check=True)
    #print("sed -i \'s/^>.*/>"+items[-1]+"/\' "+filepath2)


    subprocess.run("cp "+filepath2+" "+"${params.output}"+"/assemblies/", shell=True, check=True)   
    subprocess.run('cp ' + "${mypath}" + '/variants/' + items[-1] + '.variants.tsv ' + "${params.output}"+'/variants/', shell=True, check=True)

    filepath3="${mypath}/kraken_out/"+items[-1]+".report"
    with open(filepath3, 'r') as kreport:
        lines = kreport.readlines()
        for l in lines:
            l_parse = l.lstrip().rstrip().split("\t")
            percent = l_parse[0]
            tax_level = l_parse[3]
            tax = l_parse[5].lstrip()
            if tax_level == 'S':
                break


    with open("${mypath}"+"/report.txt", 'w') as report:
        header = ['sampleID', 'k_species', 'k_percent', 'reference', 'start', 'end', 'num_raw_reads', 'num_clean_reads', 'num_mapped_reads', 'percent_mapped_clean_reads', 'cov_bases_mapped', 'percent_genome_cov_map', 'mean_depth', 'mean_base_qual', 'mean_map_qual', 'assembly_length', 'assembly_numN', 'assembly_percent_ref_genome_cov']
        report.write('\t'.join(map(str,header)) + '\n')
        results = [items[-1], tax, percent, ref_name, start, end, raw_reads, clean_reads, reads_mapped, percent_map, cov_bases, cov, depth, baseq, mapq, num_bases, ns, pg]
        report.write('\t'.join(map(str,results)) + '\n')
    /$
}
