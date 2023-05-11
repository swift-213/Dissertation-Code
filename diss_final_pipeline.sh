#harddrive path:
cd /Volumes/Seagate/Frankie_DTOL_lep_project/

#One Drive path:
cd /Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/

#Mapping and creating the vcf plus filtered and alt versions
    #Manual
        asm=ASM10
        alt_assembly=ilOmpLuno1.1
        alt_GCA=GCA_916610225.1
        ref_assembly=ilOmpLuno1.1
        ref_GCA=GCA_916610215.1

    #Sets
        set=set6
        big_textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}_full
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            alt_GCA=$(echo ${line} | cut -d ";" -f 2 )
            ref_assembly=$(echo ${line} | cut -d ";" -f 3 )
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${alt_GCA} ${ref_assembly} ${ref_GCA}
            asm=ASM10

            #Map the alt to the reference assemblies and make bam file
            output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/minimap2_output/${alt_assembly}_alignment.bam 
            alternate=/Volumes/Seagate/Frankie_DTOL_lep_project/Fasta_files/${alt_GCA}/ncbi_dataset/data/${alt_GCA}/${alt_GCA}_${alt_assembly}_alternate_haplotype_genomic.fna
            ref=/Volumes/Seagate/Frankie_DTOL_lep_project/Fasta_files/${ref_GCA}/ncbi_dataset/data/${ref_GCA}/${ref_GCA}_${ref_assembly}_genomic.fna
            minimap2 -ax asm10 ${ref} ${alternate} | samtools view -Sb > ${output} 

            output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/minimap2_output/${alt_assembly}_alignment.bam 
            samtools_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/samtools_output/${alt_assembly}_alignment.sort.bam
            samtools sort ${output} > ${samtools_output}
            samtools index ${samtools_output}

            bcftools mpileup -f ${ref} --min-MQ 60 --annotate FORMAT/DP ${samtools_output} > ${samtools_output%.sort.bam}.mpileup
            bgzip ${samtools_output%.sort.bam}.mpileup
            bcftools index -f ${samtools_output%.sort.bam}.mpileup.gz

            #Turning VCF into bam 
            vcf_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_files/${alt_assembly}_alignment_${asm}_mcall.vcf
            bcftools call -m --ploidy 1 ${samtools_output%.sort.bam}.mpileup.gz | bgzip > ${vcf_output}
            bcftools index ${vcf_output}.gz
            

            #Filtering the VCF to DP=1 
            vcf_depth_filter=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_depth_filter/${alt_assembly}_alignment_${asm}_mcall_DP1.vcf.gz
            bcftools filter -i 'FORMAT/DP <= 1 & INFO/DP<=1' ${vcf_output} | bgzip > ${vcf_depth_filter}
            bcftools index ${vcf_depth_filter}
    
        done
        
        #Getting the numbered VCF file rather than long name Easier to do set by set due to check <- REQUIRES MANUAL INPUT TO CHECK RENAME FILE 
        set=set6
        big_textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/alex_files_loop.txt
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            alt_GCA=$(echo ${line} | cut -d ";" -f 2 )
            ref_assembly=$(echo ${line} | cut -d ";" -f 3 )
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${alt_GCA} ${ref_assembly} ${ref_GCA}

            output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/minimap2_output/${alt_assembly}_alignment.bam 
            alternate=/Volumes/Seagate/Frankie_DTOL_lep_project/Fasta_files/Fasta_files/${alt_GCA}/ncbi_dataset/data/${alt_GCA}/${alt_GCA}_${alt_assembly}_alternate_haplotype_genomic.fna
            ref=/Volumes/Seagate/Frankie_DTOL_lep_project/Fasta_files/Fasta_files/${ref_GCA}/ncbi_dataset/data/${ref_GCA}/${ref_GCA}_${ref_assembly}_genomic.fna

            vcf_depth_filter=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_depth_filter/${alt_assembly}_alignment_${asm}_mcall_DP1.vcf.gz
            vcf_depth_filter_chrom_alt=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_depth_filter/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt.vcf

            #get index, get first column
            rename=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/text_files/${alt_assembly}_rename.txt
            samtools faidx ${ref}
                
            #add row number column
            cut -f1 ${ref}.fai | awk '{print $0,NR}' | > ${rename}
            cp ${rename} /Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}/${alt_assembly}_chrom_num.txt
            
            Num_text=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}/${alt_assembly}_chrom_num.txt
            code ${Num_text}
            code ${rename}
            cat ${Num_text} | less -S
        done
    
        #use bcftools to change name in vcf file - Check rename file first!
        alt_assembly=ilOmpLuno1.1
        alt_GCA=GCA_916610225.1
        ref_assembly=ilOmpLuno1.1
        ref_GCA=GCA_916610215.1
        
        set=set6
        big_textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}_full.txt
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            alt_GCA=$(echo ${line} | cut -d ";" -f 2 )
            ref_assembly=$(echo ${line} | cut -d ";" -f 3 )
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${alt_GCA} ${ref_assembly} ${ref_GCA}
            
            vcf_depth_filter=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_depth_filter/${alt_assembly}_alignment_${asm}_mcall_DP1.vcf.gz
            vcf_depth_filter_chrom_alt=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_depth_filter/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt.vcf
            rename=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/text_files/${alt_assembly}_rename.txt
            
            bcftools annotate ${vcf_depth_filter} --rename-chrs ${rename} > ${vcf_depth_filter_chrom_alt}
            bgzip ${vcf_depth_filter_chrom_alt}
            bcftools index -f ${vcf_depth_filter_chrom_alt}.gz
        done

#Getting four fold degenerate sites:
    for i in {1..8} 
        do 
        echo $i

        number=$i
        set=set$i
        echo $set
        big_textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}_full.txt
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            alt_GCA=$(echo ${line} | cut -d ";" -f 2 )
            ref_assembly=$(echo ${line} | cut -d ";" -f 3 )
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${alt_GCA} ${ref_assembly} ${ref_GCA}
            asm=ASM10
        
            code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/genomics_general/
            rename=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/text_files/${alt_assembly}_rename.txt
            rename_four=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/text_files/${alt_assembly}_rename.txt
            gff=/Volumes/Seagate/Frankie_DTOL_lep_project/gff/files/${ref_GCA}-genes.gff3.gz
            ref=/Volumes/Seagate/Frankie_DTOL_lep_project/Fasta_files/${ref_GCA}/ncbi_dataset/data/${ref_GCA}/${ref_GCA}_${ref_assembly}_genomic.fna
            vcf_depth_filter_chrom_alt=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_depth_filter/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt.vcf.gz
            four_D_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/fourD_sites/${alt_assembly}_coding_site_types.tsv.gz
            four_D_bed_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/fourD_sites/${alt_assembly}_fourDsites.bed

            python3 ${code_path}codingSiteTypes.py -a ${gff} -f gff3 -r ${ref} -v ${vcf_depth_filter_chrom_alt} --scaffoldLookup ${rename_four} --ignoreConflicts --useAnnotationScaffoldNames | bgzip > ${four_D_output}
            gunzip -c  ${four_D_output} | awk 'BEGIN {OFS="\t"}; $5=="4" {print($1"\t"$2-1"\t"$2)}'|  > ${four_D_bed_output}

            output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/fourD_sites/vcfs/${alt_assembly}_fourDsites.vcf
            bedtools intersect -wa -a ${vcf_depth_filter_chrom_alt} -b ${four_D_bed_output} -header > $output 
        done

        #getting 4D snp positions:
        big_textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/alex_files_loop.txt
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            set=$(echo ${line} | cut -d ";" -f 2 )
            echo ${alt_assembly} ${set}
            asm=ASM10
            code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project
            vcf_fourD=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/fourD_sites/vcfs/${alt_assembly}_fourDsites.vcf
            Assembly_info=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/fourD_sites/snps/${alt_assembly}_fourD_snps.txt
            chrom_text=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${alt_assembly}_chrom_num.txt
        
            
            python3 ${code_path}/fourD_snp_finder.py -i ${vcf_fourD} -i2 ${chrom_text} -o ${Assembly_info} 


            Assembly_info=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/fourD_sites/snps/${alt_assembly}_fourD_snps.txt
            Assembly_info_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/fourD_sites/snps/${alt_assembly}_fourD_snps_sed.txt
            sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}


#Whole assembly analysis
    for i in {1..8} 
    do 
    echo $i

    number=$i
    set=set$i
    echo $set
    set=set5
    textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}.txt
    asm=ASM10
    for line in $(cat ${textfile})
        do 
        alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
        echo ${alt_assembly}
   
        code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project
        vcf_depth_filter_chrom_alt=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_depth_filter/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt.vcf.gz
        Assembly_info=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/Assembly_info/${alt_assembly}_assembly_info.txt
        indel_finder=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelFinder/${alt_assembly}_IndelFinder_no_CDS_trunc_new.txt
        indel_binner_100=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelBinner/hundred/${alt_assembly}_IndelBinner_100_no_CDS_trunc.txt
        indel_binner_10000=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelBinner/ten_thousand/${alt_assembly}_IndelBinner_10000_no_CDS_trunc.txt
        chrom_text=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}/${alt_assembly}_chrom_num.txt
        every_indel_length=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelBinner/all_indel_lengths/${alt_assembly}_IndelBinner_all_pos.txt
        
        python3 ${code_path}/full_script_reborn.py -i ${vcf_depth_filter_chrom_alt} -i2 ${chrom_text} -o ${Assembly_info} -o2 ${indel_finder} -o3 ${indel_binner_100} -o4 ${indel_binner_10000} -o5 ${every_indel_length}


        Assembly_info=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/Assembly_info/${alt_assembly}_assembly_info.txt
        Assembly_info_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/Assembly_info/${alt_assembly}_assembly_info_sed.txt
        sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}

        indel_finder=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelFinder/${alt_assembly}_IndelFinder_new.txt
        indel_finder_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelFinder/${alt_assembly}_IndelFinder_new_sed.txt
        sed 's/\[//g' ${indel_finder} | sed 's/\]//g' | > ${indel_finder_sed}

      

        indel_binner_100=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelBinner/hundred/${alt_assembly}_IndelBinner_100_no_CDS_trunc.txt
        indel_binner_100_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelBinner/hundred/${alt_assembly}_IndelBinner_100_sed.txt
        sed 's/{//g' ${indel_binner_100} | sed 's/}//g' | sed 's/:/,/g' | sed 's/'\''//g' | > ${indel_binner_100_sed}
        
        indel_binner_10000=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelBinner/ten_thousand/${alt_assembly}_IndelBinner_10000_no_CDS_trunc.txt
        indel_binner_10000_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelBinner/ten_thousand/${alt_assembly}_IndelBinner_10000_no_CDS_trunc_sed.txt
        sed 's/{//g' ${indel_binner_10000} | sed 's/}//g' | sed 's/:/,/g' | sed 's/'\''//g' | > ${indel_binner_10000_sed}

        every_indel_length=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelBinner/all_indel_lengths/${alt_assembly}_IndelBinner_all_pos.txt
        every_indel_length_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelBinner/all_indel_lengths/${alt_assembly}_IndelBinner_all_pos_sed.txt
        sed 's/{//g' ${every_indel_length} | sed 's/}//g' | sed 's/:/,/g' | sed 's/'\''//g' | > ${every_indel_length_sed}
    done
    done
#CDS analysis
    #Pull CDS pos out of gff file:
        for i in {1..8} 
        do 
        echo $i

        number=$i
        set=set$i
        echo $set
        big_textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}_full.txt

        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${ref_GCA}
       
            gff=/Volumes/Seagate/Frankie_DTOL_lep_project/gff/files/${ref_GCA}-genes.gff3
            final_file=/Volumes/Seagate/Frankie_DTOL_lep_project/gff/exon_only_files/${alt_assembly}_CDS_info.bed
            final_file_short=/Volumes/Seagate/Frankie_DTOL_lep_project/gff/exon_only_files/${alt_assembly}_CDS.bed
            final_file_sorted=/Volumes/Seagate/Frankie_DTOL_lep_project/gff/exon_only_files/${alt_assembly}_CDS_info_sorted.bed
            final_file_sorted_merged=/Volumes/Seagate/Frankie_DTOL_lep_project/gff/exon_only_files/${alt_assembly}_CDS_info_sorted_merged.bed
            final_file_sorted_short=/Volumes/Seagate/Frankie_DTOL_lep_project/gff/exon_only_files/${alt_assembly}_CDS_info_sorted_short.bed
            final_file_sorted_merged_short=/Volumes/Seagate/Frankie_DTOL_lep_project/gff/exon_only_files/${alt_assembly}_CDS_info_sorted_merged_short.bed
            #pull out all the exon rows
                awk '{if($3=="CDS") print $0}' ${gff} > all_exons.bed
            #remove and append column 1 and 9 - look at rishi code and change the chrom names in one of the files so they match 
                awk '{if($1!="Z") print $0}' all_exons.bed > all_exons_1.bed
                awk '{if($1!="W") print $0}' all_exons_1.bed > all_exons_2.bed
                awk '{if($1!="MT") print $0}' all_exons_2.bed > all_exons_3.bed
                awk '{print $1}' all_exons_3.bed  | > our_column1.bed

            #Getting only geneID from column 9 
                awk '{print $9}' all_exons_3.bed | cut -d ';' -f 1 | > our_column9.bed
                awk '{print $4-1}' all_exons_3.bed > all_exons_start.bed
                awk '{print $5}' all_exons_3.bed > all_exons_end.bed

            #Add it all together then sort and merge so the exons are whole
                paste our_column1.bed all_exons_start.bed all_exons_end.bed our_column9.bed > ${final_file}
                paste our_column1.bed all_exons_start.bed all_exons_end.bed > ${final_file_short}
                sort -k1,1 -k2,2n ${final_file} > ${final_file_sorted}
                bedtools merge -i ${final_file_sorted} -c 4 -o distinct > ${final_file_sorted_merged}
                sort -k1,1 -k2,2n ${final_file_short} > ${final_file_sorted_short}
                bedtools merge -i ${final_file_sorted_short} > ${final_file_sorted_merged_short}

    #SNP CDS analysis
        #Making the CDS included and excluded files for the vcf with the number chrom 
            textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/all_gff_sets.txt
            for line in $(cat ${textfile})
                do 
                alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                set=$(echo ${line} | cut -d ";" -f 2 )
                echo ${alt_assembly} ${set}

            
                asm=ASM10
                exon_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/gff/exon_only_files/${alt_assembly}_CDS_info_sorted_merged_short.bed
                vcf_depth_filter_chrom_alt=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_depth_filter/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt.vcf.gz
                output_include=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt_CDS_in.vcf
                output_exclude=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt_CDS_out.vcf

                bcftools view ${vcf_depth_filter_chrom_alt} --targets-file ${exon_file} -o ${output_include}
                bgzip -f ${output_include}
                bcftools index ${output_include}.gz
                
                #Getting snps for whole cds regions - need for indel diversity denominator 
                code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project
                output_include=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt_CDS_in.vcf
                Assembly_info=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_cds_total_snps.txt
                chrom_text=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}/${alt_assembly}_chrom_num.txt
            
                
                python3 ${code_path}/fourD_snp_finder.py -i ${output_include}.gz -i2 ${chrom_text} -o ${Assembly_info} 

                Assembly_info=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_cds_total_snps.txt
                Assembly_info_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_cds_total_snps_sed.txt
                sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}


                #getting four d positions and then snps at the four d positions 
                four_D_bed_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/fourD_sites/${alt_assembly}_fourDsites.bed
                output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_fourD_sites_CDS.vcf
                bedtools intersect -wa -a ${output_include}.gz -b ${four_D_bed_output} -header > $output 
                

                code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project
                vcf_fourD=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_fourD_sites_CDS.vcf
                Assembly_info=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_fourD_snps.txt
                chrom_text=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}/${alt_assembly}_chrom_num.txt
            
                
                python3 ${code_path}/fourD_snp_finder.py -i ${vcf_fourD} -i2 ${chrom_text} -o ${Assembly_info} 


                Assembly_info=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_fourD_snps.txt
                Assembly_info_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_fourD_snps_sed.txt
                sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}


                #making vcfs from the repeat regions
                Repeat_positions_vcf=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/Repeat_regions/final/${alt_assembly}_repeat_regions_sorted.bed
                vcf_depth_filter_chrom_alt=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_depth_filter/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt.vcf.gz
                repeats_output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/Repeat_regions/VCF/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt_repeats.vcf
                
                bcftools view ${vcf_depth_filter_chrom_alt} --targets-file ${Repeat_positions_vcf} -o ${repeats_output}
                bgzip -f ${repeats_output}
                bcftools index ${repeats_output}.gz


                #Getting snps for whole repeat regions - need for indel diversity denominator 
                code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project
                repeats_output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/Repeat_regions/VCF/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt_repeats.vcf
                Assembly_info=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/Repeat_regions/${alt_assembly}_repeats_total_snps.txt
                chrom_text=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}/${alt_assembly}_chrom_num.txt
            
                
                python3 ${code_path}/fourD_snp_finder.py -i ${repeats_output}.gz -i2 ${chrom_text} -o ${Assembly_info} 

                Assembly_info=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/Repeat_regions/${alt_assembly}_repeats_total_snps.txt
                Assembly_info_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/Repeat_regions/${alt_assembly}_repeats_total_snps_sed.txt
                sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}



                four_D_bed_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/fourD_sites/${alt_assembly}_fourDsites.bed
                output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/Repeat_regions/VCF/${alt_assembly}_fourD_sites_CDS.vcf
                bedtools intersect -wa -a ${repeats_output}.gz -b ${four_D_bed_output} -header > $output 

                code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project
                vcf_fourD=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/Repeat_regions/VCF/${alt_assembly}_fourD_sites_CDS.vcf
                Assembly_info=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/Repeat_regions/VCF/${alt_assembly}_fourD_snps.txt
                chrom_text=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}/${alt_assembly}_chrom_num.txt
            
                
                python3 ${code_path}/fourD_snp_finder.py -i ${vcf_fourD} -i2 ${chrom_text} -o ${Assembly_info} 


                Assembly_info=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/Repeat_regions/VCF/${alt_assembly}_fourD_snps.txt
                Assembly_info_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/Repeat_regions/VCF/${alt_assembly}_fourD_snps_sed.txt
                sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}


                #intergenic get estimates in R by subtracting the other values from our total ones
        
        #pulling out data for scatterplots of SNP and Indel diveristy - need a gff
            set=set3
            textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}.txt
            for line in $(cat ${textfile})
                do 
                alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                echo ${alt_assembly}
            
                asm=ASM10
                exon_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/gff/exon_only_files/${alt_assembly}_exon_info_sorted_merged.bed
                vcf_depth_filter_chrom_alt=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/vcf_depth_filter/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt.vcf.gz
                output_include=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt_CDS_in.vcf.gz
                output_exclude=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt_CDS_out.vcf.gz


                code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/
                output_include=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt_CDS_in.vcf.gz
                Assembly_info=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_included_SNP/Assembly_info/${alt_assembly}_assembly_info.txt
                indel_finder=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelFinder/${alt_assembly}_IndelFinder_CDS_trunc.txt
                indel_binner_100=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelBinner/hundred/${alt_assembly}_IndelBinner_100_CDS_trunc.txt
                indel_binner_10000=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelBinner/ten_thousand/${alt_assembly}_IndelBinner_10000_CDS_trunc.txt
                every_indel_length=
                chrom_text=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}/${alt_assembly}_chrom_num.txt

                python3 ${code_path}full_script_reborn.py -i ${output_include} -i2 ${chrom_text} -o ${Assembly_info} -o2 ${indel_finder} -o3 ${indel_binner_100} -o4 ${Indel_binner_10000}


                code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/
                output_exclude=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/${alt_assembly}_alignment_${asm}_mcall_DP1_chrom_alt_CDS_out.vcf.gz
                Assembly_info=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/Assembly_info/${alt_assembly}_assembly_info.txt
                indel_finder=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelFinder/${alt_assembly}_IndelFinder_no_CDS_trunc.txt
                indel_binner_100=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelBinner/hundred/${alt_assembly}_IndelBinner_100_no_CDS_trunc.txt
                indel_binner_10000=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelBinner/ten_thousand/${alt_assembly}_IndelBinner_10000_no_CDS_trunc.txt
                chrom_text=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}/${alt_assembly}_chrom_num.txt

                python3 ${code_path}/full_script_reborn.py -i ${output_exclude} -i2 ${chrom_text} -o ${Assembly_info} -o2 ${indel_finder} -o3 ${indel_binner_100} -o4 ${Indel_binner_10000}
            #Sorting of the outputs for R 
                #Exons Included:
                    asm=ASM10
                    #sorting indel pos output - can then be loaded into R as a table before R then calculates the length of them! 
                        sed_indel_pos_input=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelPos/${alt_assembly}_Indelpos_CDS_trunc.txt
                        indel_pos_sed_output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelPos/${alt_assembly}_Indelpos_CDS_sed_trunc.txt
                        sed 's/\[//g' ${sed_indel_pos_input} |  sed 's/\]//g' |sed 's/\,//g' | sed "s/'//g" | > ${indel_pos_sed_output}
                    
                    #Sorting indel finder output - just removes brackets and stuff for loading into R
                        input_sed_indelfinder=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelFinder/${alt_assembly}_IndelFinder_CDS_trunc.txt 
                        IndelFinder_sed_output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelFinder/${alt_assembly}_IndelFinder_CDS_sed_trunc.txt
                        sed 's/\[//g' ${input_sed_indelfinder} | sed 's/\]//g' | head > ${IndelFinder_sed_output}
                        
                    #Sorting <100 indel binning output - removes brackets, gets column 1 to all the same decimal places so they can be sorted correctly and removes empty columns - also removes the > and < columns but these can still be seen in the orginal output
                        indel_binner_output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelBinner/100/${alt_assembly}_IndelBinner_100_CDS_trunc.txt
                        indel_binner_output_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelBinner/100/${alt_assembly}_IndelBinner_100_CDS_sed_trunc.txt
                        indel_binner_output_sed_r=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelBinner/100/${alt_assembly}_IndelBinner_100_CDS_sed_r_trunc.txt

                        sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > save.txt
                        awk '{print $1}' save.txt |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                        awk '{print $2}' save.txt | > col_2.txt
                        paste col_1.txt col_2.txt > unsorted.txt
                        awk '$1!=0000' unsorted.txt > sorted.txt
                        sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                        sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r} 

                    #sorting indel binner up to 10,000 - more complex than <100 output as column 1 needs more work but does the same thing pretty much - also removes the > and < columns but these can still be seen in the orginal output 
                        indel_binner_output_10000=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelBinner/10000/${alt_assembly}_IndelBinner_10000_CDS_trunc.txt
                        indel_binner_output_10000_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_included_SNP/IndelBinner/10000/${alt_assembly}_IndelBinner_10000_exons_CDS_trunc.txt
                        sed 's/{//g' ${indel_binner_output_10000} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > save.txt
                        awk '{print $1}' save.txt |  sed 's/-/'\ '/g' |  awk '{printf ("%05d\n", $1)}' | > column1_1.txt
                        awk '{print $1}' save.txt |  sed 's/-/'\ '/g' |  awk '{printf ("%05d\n", $2)}' | > column1_2.txt
                        paste column1_1.txt column1_2.txt > column1_3.txt
                        awk ' { print $0, $1 "-" $NF } ' column1_3.txt | awk '{print $3}' | > column1.txt
                        awk '{print $2}' save.txt | > column2.txt
                        paste column1.txt column2.txt > unsorted.txt
                        awk '$1!=00000-00000' unsorted.txt > sorted.txt
                        sort -k 1 sorted.txt > ${indel_binner_output_10000_sed} 

                

                #Exons not included
                    sed_indel_pos_input=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelPos/${alt_assembly}_Indelpos_no_CDS_trunc.txt
                    indel_pos_sed_output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelPos/${alt_assembly}_Indelpos_no_CDS_sed_trunc.txt
                    sed 's/\[//g' ${sed_indel_pos_input} |  sed 's/\]//g' |sed 's/\,//g' | sed "s/'//g" | > ${indel_pos_sed_output}

                    #Sorting indel finder output - just removes brackets and stuff for loading into R
                        input_sed_indelfinder=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelFinder/${alt_assembly}_IndelFinder_no_CDS_trunc.txt 
                        IndelFinder_sed_output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelFinder/${alt_assembly}_IndelFinder_no_CDS_sed_trunc.txt

                        sed 's/\[//g' ${input_sed_indelfinder} | sed 's/\]//g' | head > ${IndelFinder_sed_output}

                    #Sorting <100 indel binning output - removes brackets, gets column 1 to all the same decimal places so they can be sorted correctly and removes empty columns - also removes the > and < columns but these can still be seen in the orginal output
                        indel_binner_output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelBinner/100/${alt_assembly}_IndelBinner_100_no_CDS_trunc.txt
                        indel_binner_output_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelBinner/100/${alt_assembly}_IndelBinner_100_no_CDS_sed_trunc.txt
                        indel_binner_output_sed_r=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelBinner/100/${alt_assembly}_IndelBinner_100_no_CDS_sed_r_trunc.txt

                        sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > save.txt
                        awk '{print $1}' save.txt |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                        awk '{print $2}' save.txt | > col_2.txt
                        paste col_1.txt col_2.txt > unsorted.txt
                        awk '$1!=0000' unsorted.txt > sorted.txt
                        sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                        sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r} 

                    #sorting indel binner up to 10,000 - more complex than <100 output as column 1 needs more work but does the same thing pretty much - also removes the > and < columns but these can still be seen in the orginal output
                        indel_binner_output_10000=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelBinner/10000/${alt_assembly}_IndelBinner_10000_no_CDS_trunc.txt
                        indel_binner_output_10000_sed=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/IndelBinner/10000/${alt_assembly}_IndelBinner_10000_no_CDS_sed_trunc.txt

                        sed 's/{//g' ${indel_binner_output_10000} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > save.txt
                        awk '{print $1}' save.txt |  sed 's/-/'\ '/g' |  awk '{printf ("%05d\n", $1)}' | > column1_1.txt
                        awk '{print $1}' save.txt |  sed 's/-/'\ '/g' |  awk '{printf ("%05d\n", $2)}' | > column1_2.txt
                        paste column1_1.txt column1_2.txt > column1_3.txt
                        awk ' { print $0, $1 "-" $NF } ' column1_3.txt | awk '{print $3}' | > column1.txt
                        awk '{print $2}' save.txt | > column2.txt
                        paste column1.txt column2.txt > unsorted.txt
                        awk '$1!=00000-00000' unsorted.txt > sorted.txt
                        sort -k 1 sorted.txt > ${indel_binner_output_10000_sed} 

    #Indel analysis
        #Intersect to get indels in cds
            for i in {1..8} 
            do 
            echo $i

            number=$i
            set=set$i
            echo $set
            big_textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}.txt

            for line in $(cat ${big_textfile})
                do 
                alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                echo ${alt_assembly}

                #Getting bed files rather than text
                    indel_pos_sed_output_text=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/outputs/IndelFinder/${alt_assembly}_IndelFinder_new_sed.txt
                    rename_text_sed2=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/IndelFinder/${alt_assembly}_IndelFinder_new_sed.bed
                    cp ${indel_pos_sed_output_text} ${rename_text_sed2}
                
                    exon_file_short=/Volumes/Seagate/Frankie_DTOL_lep_project/gff/exon_only_files/${alt_assembly}_CDS_info_sorted_merged_short.bed
                    exon_file_short_sed=/Volumes/Seagate/Frankie_DTOL_lep_project/gff/exon_only_files/${alt_assembly}_CDS_info_sorted_merged_short_sed.bed
                    indel_pos_sed_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/IndelFinder/${alt_assembly}_IndelFinder_new_sed.bed
                    indel_pos_sed_output_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/IndelFinder/${alt_assembly}_IndelFinder_new_sed_ready.bed
                    output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/IndelOverlapExon/raw/${alt_assembly}_indels_intersect_CDS_neww.bed
                
                #Bedtools hates spaces so have to makesure input is tabbed and also that the indel pos jave 0 based column 2 and 1 based column 3:
                    sed 's/ /\t/g' ${exon_file_short} > ${exon_file_short_sed}
                    
                    awk '{print $2-1}' $rename_text_sed2 > column_2.txt
                    awk '{print $1}' $rename_text_sed2 > column_1.txt
                    awk '{print $3}' $rename_text_sed2 > column_3.txt
                    
                    paste column_1.txt column_2.txt column_3.txt > ${indel_pos_sed_output}
                    sed 's/ /\t/g' ${indel_pos_sed_output} | sed 's/,//g' |sed 's/'\''//g' | > ${indel_pos_sed_output_sed}

                #Actual intersect command 
                    bedtools intersect -wao -a ${indel_pos_sed_output_sed} -b ${exon_file_short_sed} > ${output}

        #R -STEPS
            #this gets the assemblies ready to go through the indel dictionary sorter from r
            #Indels in CDS:
                text_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/all_gff_shell.txt
                for line in $(cat ${text_file})
                    do 
                    alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                    echo ${alt_assembly} 

                    indels_in_CDS_length=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_in_exons/${alt_assembly}_Indels_in_CDS_trunc_raw_needs_sorted_new.csv
                    output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_in_exons/${alt_assembly}_Indels_in_CDS_trunc_raw_needs_sorted_sed_new.csv

                    awk 'NR!=1{print $1}' ${indels_in_CDS_length} > ${output}

                #Indel dictionary sorter from R
                    text_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/all_gff_shell.txt
                    for line in $(cat ${text_file})
                        do 
                        alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                        echo ${alt_assembly} 
                        code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/
                        input=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_in_exons/${alt_assembly}_Indels_in_CDS_trunc_raw_needs_sorted_sed_new.csv
                        dict_100=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_included_SNP/R_indel_100_binner/${alt_assembly}_100_binner_r_indels_in_CDS_new.txt
                        greater_than_100_dict=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_included_SNP/R_indel_greater_than_100_binner/${alt_assembly}_greater_than_100_binner_r_indels_in_CDS_new.txt

                        python3 ${code_path}/exon_R_indel_sorter.py -i ${input} -o ${dict_100} -o2 ${greater_than_100_dict}

                        indel_binner_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_included_SNP/R_indel_100_binner/${alt_assembly}_100_binner_r_indels_in_CDS_new.txt
                        indel_binner_output_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_included_SNP/R_indel_100_binner/${alt_assembly}_100_binner_r_indels_in_CDS_sed_new.txt
                        indel_binner_output_sed_r=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_included_SNP/R_indel_100_binner/${alt_assembly}_100_binner_r_indels_in_CDS_sed_r_new.txt

                        sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > save.txt
                        awk '{print $1}' save.txt |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                        awk '{print $2}' save.txt | > col_2.txt
                        paste col_1.txt col_2.txt > unsorted.txt
                        awk '$1!=0000' unsorted.txt > sorted.txt
                        sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                        sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r} 

                #Indels outside CDS - getting lengths and disctionary sorting is all in one code file!
                text_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/all_gff_shell.txt
                    for line in $(cat ${text_file})
                        do 
                        alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                        echo ${alt_assembly} 
                
                input_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_not_in_exons/${alt_assembly}_Indels_not_in_CDS_trunc_raw_needs_sorted_new.csv
                input_file_2=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_not_in_exons/${alt_assembly}_Indels_not_in_CDS_trunc_raw_sorted_new.csv
                sed 's/" //g' $input_file | sed 's/"//g' | sed 's/ //g' | sed 's/\t/,/g' | > $input_file_2


                    text_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/all_gff_shell.txt
                    for line in $(cat ${text_file})
                        do 
                        alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                        echo ${alt_assembly} 
                        code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/
                        input_file_2=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_not_in_exons/${alt_assembly}_Indels_not_in_CDS_trunc_raw_sorted_new.csv
                        indel_total_lengths=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_not_in_exons/${alt_assembly}_Indels_not_in_CDS_trunc_raw_sorted_lengths_new.csv
                        dict_100=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/R_indel_100_binner/${alt_assembly}_100_binner_r_indels_not_in_CDS_new.txt
                        greater_than_100_dict=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/R_indel_greater_than_100_binner/${alt_assembly}_greater_than_100_binner_r_indels_not_in_CDS_new.txt
                        Number_of_indels=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_not_in_exons/${alt_assembly}_number_of_indels_new.txt

                        python3 ${code_path}/Indels_outside_CDS_counter.py -i ${input_file_2} -o ${indel_total_lengths} -o2 ${dict_100} -o3 ${greater_than_100_dict} -o4 ${Number_of_indels}

                        indel_binner_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/R_indel_100_binner/${alt_assembly}_100_binner_r_indels_not_in_CDS_new.txt
                        indel_binner_output_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/R_indel_100_binner/${alt_assembly}_100_binner_r_indels_not_in_CDS_sed_new.txt
                        indel_binner_output_sed_r=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/exons_removed_SNP/R_indel_100_binner/${alt_assembly}_100_binner_r_indels_not_in_CDS_sed_r_new.txt
                        sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > save.txt
                        awk '{print $1}' save.txt |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                        awk '{print $2}' save.txt | > col_2.txt
                        paste col_1.txt col_2.txt > unsorted.txt
                        awk '$1!=0000' unsorted.txt > sorted.txt
                        sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                        sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r}
                        
                        Number_of_indels=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_not_in_exons/${alt_assembly}_number_of_indels_new.txt
                        Number_of_indels_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_not_in_exons/${alt_assembly}_number_of_indels_sed_new.txt
                        sed 's/{//g' ${Number_of_indels} | sed 's/}//g' | sed 's/'\''//g' | sed 's/Number_of_indels://g' | > ${Number_of_indels_sed}




#Repeats analysis 
    #Repeat_Regions, runs through the reference fasta and pulls out the positions of repeat regions
            for i in {7..8} 
            do 
            echo $i

            number=$i
            set=set$i
            echo $set
        
        big_textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}_full.txt
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            ref_assembly=$(echo ${line} | cut -d ";" -f 3)
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${ref_assembly} ${ref_GCA}



            code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/
            chrom_names_R=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/${set}/${alt_assembly}_chrom_R.txt
            ref=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/Fasta_files/${ref_GCA}/ncbi_dataset/data/${ref_GCA}/${ref_GCA}_${ref_assembly}_genomic.fna
            Repeat_positions=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/Repeat_regions/${alt_assembly}_repeat_regions.bed
            Repeat_positions_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/Repeat_regions/${alt_assembly}_repeat_regions_sed.bed
            #repeat regions come out with correct bed file formatiing (col2=0 col3=1 based)

            python3 ${code_path}/Repeat_finder.py -i ${ref} -i2 ${chrom_names_R} -o ${Repeat_positions} 
            sed 's/\[//g' ${Repeat_positions} | sed 's/\]//g' | > ${Repeat_positions_sed}
    cat $Repeat_positions_sed | less -S
    #Intersecting the repeat regions with everything else regions 
    #Edit bed files so that column 1 = 0 based and column 2 = 1 based (-1 from all first columns), awk pull out $2 and then minus one from it in order to get the 0 based positions.
     for i in {1..8} 
            do 
            echo $i

            number=$i
            set=set$i
            echo $set
     #set=set8
        big_textfile=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/all_gff_shell.txt
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            echo ${alt_assembly} 

            #You dont need to do $2 minus one here as the indel pos comes from CDS analysis where they are already in bed file form!
            
            indels_not_in_CDS=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_not_in_exons/${alt_assembly}_Indels_not_in_CDS_trunc_raw_needs_sorted_new.csv
            
            
            indels_not_in_CDS_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/Exons/Indels_not_in_exons/${alt_assembly}_Indels_not_in_CDS_trunc_raw_needs_sorted_new_sed.bed
            awk 'NR!=1{print $0}' ${indels_not_in_CDS} | sed 's/ //g' | sed 's/"//g' | > ${indels_not_in_CDS_sed}
            

          
            Repeat_positions_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/Repeat_regions/${alt_assembly}_repeat_regions_sed.bed
            repeat_pos_ready=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/Repeat_regions/${alt_assembly}_repeat_regions_sed_ready.bed
            seddy_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/Repeat_regions/${alt_assembly}_repeats_alt_Chrom.bed
            repeat_pos_final=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/Repeat_regions/final/${alt_assembly}_repeat_regions_sorted.bed
            rename=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/text_files/${alt_assembly}_rename.txt
            
            sed 's/'\''//g' $Repeat_positions_sed | sed 's/,//g' | > ${repeat_pos_ready}

            join ${repeat_pos_ready} ${rename} > ${seddy_output}
            awk '{print $4"\t"$2"\t"$3}' $seddy_output | > ${repeat_pos_final}
           
    
            output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/IndelOverlapRepeat/${alt_assembly}_repeats_intersect_indels.bed
            bedtools intersect -wao -a ${indels_not_in_CDS_sed} -b ${repeat_pos_final} > ${output}

            
    #just use indels with no repeats in them at all?

        #Python analysis - in repeats - go to R and get the inputs first 
        text_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/all_gff_shell.txt
                for line in $(cat ${text_file})
                    do 
                    alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                    echo ${alt_assembly} 

                    input=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/R_results/Repeat_indels/Indels_fully_in_repeats/${alt_assembly}_Indels_in_repeats_trunc_raw_needs_sorted_new.csv
                    output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/R_results/Repeat_indels/Indels_fully_in_repeats/${alt_assembly}_Indels_in_repeats_trunc_raw_sorted_new.csv

                    awk 'NR!=1{print $0}' ${input} | sed 's/ //g' | sed 's/""/"/g' | sed 's/"/,/g' | sed 's/,$//g' | > ${output}
                    
                text_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/all_gff_shell.txt
                for line in $(cat ${text_file})
                    do 
                    alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                    echo ${alt_assembly} 

                    code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/
                    input_file=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/R_results/Repeat_indels/Indels_fully_in_repeats/${alt_assembly}_Indels_in_repeats_trunc_raw_sorted_new.csv
                    indel_total_lengths=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_in_repeats/${alt_assembly}_Indels_in_repeats_trunc_raw_needs_sorted_lengths_new.csv
                    dict_100=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_in_repeats/indel_100_binner/${alt_assembly}_100_binner_r_indels_in_repeats_new.txt
                    greater_than_100_dict=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_in_repeats/Indels_greater_than_100_binner/${alt_assembly}_greater_than_100_binner_r_indels_in_repeats_new.txt
                    Number_of_indels=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_in_repeats/${alt_assembly}_number_of_indels_new.txt

                    python3 ${code_path}/Indels_outside_CDS_counter.py -i ${input_file_2} -o ${indel_total_lengths} -o2 ${dict_100} -o3 ${greater_than_100_dict} -o4 ${Number_of_indels}

                    indel_binner_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_in_repeats/indel_100_binner/${alt_assembly}_100_binner_r_indels_in_repeats_new.txt
                    indel_binner_output_sed1=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_in_repeats/indel_100_binner/${alt_assembly}_100_binner_r_indels_in_repeats_new_sed1.txt
                    indel_binner_output_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_in_repeats/indel_100_binner/${alt_assembly}_100_binner_r_indels_in_repeats_new_sed.txt
                    indel_binner_output_sed_r=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_in_repeats/indel_100_binner/${alt_assembly}_100_binner_r_indels_in_repeats_new_sed_r.txt
                    sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > ${indel_binner_output_sed1}
                    awk '{print $1}' ${indel_binner_output_sed1} |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                    awk '{print $2}' ${indel_binner_output_sed1} | > col_2.txt
                    paste col_1.txt col_2.txt > unsorted.txt
                    awk '$1!=0000' unsorted.txt > sorted.txt
                    sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                    sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r}
                    
                    cat $indel_binner_output_sed_r | less -S


                    Number_of_indels=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_in_repeats/${alt_assembly}_number_of_indels_new.txt
                    Number_of_indels_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_in_repeats/${alt_assembly}_number_of_indels_new_sed.txt
                    sed 's/{//g' ${Number_of_indels} | sed 's/}//g' | sed 's/'\''//g' | sed 's/Number_of_indels://g' | > ${Number_of_indels_sed}

                    indel_total_lengths=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_in_repeats/${alt_assembly}_Indels_in_repeats_trunc_raw_needs_sorted_lengths_new.csv
                    indel_total_lengths_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_in_repeats/${alt_assembly}_Indels_in_repeats_trunc_raw_needs_sorted_lengths_new_sed.csv
                    sed 's/\[//g' $indel_total_lengths | sed 's/\]//g' | sed 's/'\''//g' | >  $indel_total_lengths_sed 


       
        awk '{if ($3-$2 == 1) print $0}' $indels_not_in_CDS_sed | head
        #Python analysis - out repeats

                    text_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/all_gff_shell.txt
                for line in $(cat ${text_file})
                    do 
                    alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                    echo ${alt_assembly} 

                    input=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/R_results/Repeat_indels/Indels_not_in_repeats/${alt_assembly}_Indels_not_in_repeats_trunc_raw_needs_sorted_new.csv
                    output=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/R_results/Repeat_indels/Indels_not_in_repeats/${alt_assembly}_Indels_not_in_repeats_trunc_raw_sorted_new.csv


                    awk 'NR!=1{print $0}' ${input} | sed 's/ //g' | sed 's/""/"/g' | sed 's/"/,/g' | sed 's/,$//g' | > ${output}
                    
                    
                    text_file=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/text_files/all_gff_shell.txt
                    for line in $(cat ${text_file})
                        do 
                        alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                        echo ${alt_assembly} 
                        code_path=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/
                        input_file=/Users/frankieswift/Library/CloudStorage/OneDrive-Personal/Uni/4th_year/Honours_project/R_results/Repeat_indels/Indels_not_in_repeats/${alt_assembly}_Indels_not_in_repeats_trunc_raw_sorted_new.csv
                        indel_total_lengths=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_not_in_repeats/${alt_assembly}_Indels_in_repeats_trunc_raw_needs_sorted_lengths_new.csv
                        dict_100=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_not_in_repeats/indel_100_binner/${alt_assembly}_100_binner_r_indels_not_in_repeats_new.txt
                        greater_than_100_dict=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_not_in_repeats/Indels_greater_than_100_binner/${alt_assembly}_greater_than_100_binner_r_indels_not_in_repeats_new.txt
                        Number_of_indels=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_not_in_repeats/${alt_assembly}_number_of_indels_new.txt

                        python3 ${code_path}/Indels_outside_CDS_counter.py -i ${input_file} -o ${indel_total_lengths} -o2 ${dict_100} -o3 ${greater_than_100_dict} -o4 ${Number_of_indels}

                        indel_binner_output=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_not_in_repeats/indel_100_binner/${alt_assembly}_100_binner_r_indels_not_in_repeats_new.txt
                        indel_binner_output_sed1=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_not_in_repeats/indel_100_binner/${alt_assembly}_100_binner_r_indels_not_in_repeats_new_sed1.txt
                        indel_binner_output_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_not_in_repeats/indel_100_binner/${alt_assembly}_100_binner_r_indels_not_in_repeats_new_sed.txt
                        indel_binner_output_sed_r=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/outputs/repeats/indels_not_in_repeats/indel_100_binner/${alt_assembly}_100_binner_r_indels_not_in_repeats_new_r.txt
                        sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > ${indel_binner_output_sed1}
                        awk '{print $1}' ${indel_binner_output_sed1} |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                        awk '{print $2}' ${indel_binner_output_sed1} | > col_2.txt
                        paste col_1.txt col_2.txt > unsorted.txt
                        awk '$1!=0000' unsorted.txt > sorted.txt
                        sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                        sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r}
                        
                        Number_of_indels=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_not_in_repeats/${alt_assembly}_number_of_indels_new.txt
                        Number_of_indels_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_not_in_repeats/${alt_assembly}_number_of_indels_new_sed.txt
                        sed 's/{//g' ${Number_of_indels} | sed 's/}//g' | sed 's/'\''//g' | sed 's/Number_of_indels://g' | > ${Number_of_indels_sed}
        

                        indel_total_lengths=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_not_in_repeats/${alt_assembly}_Indels_in_repeats_trunc_raw_needs_sorted_lengths_new.csv
                        indel_total_lengths_sed=/Users/frankieswift/OneDrive/Uni/4th_year/Honours_project/R_results/repeats/Indels_not_in_repeats/${alt_assembly}_Indels_in_repeats_trunc_raw_needs_sorted_lengths_new_sed.csv
                        sed 's/\[//g' $indel_total_lengths | sed 's/\]//g' | sed 's/'\''//g' | >  $indel_total_lengths_sed 