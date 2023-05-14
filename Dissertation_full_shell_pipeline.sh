#A full automated pipeline of the project. File paths are arbritaray and would differ depending on the location of the files on the computer the code was run on. This is a blank version to show the pipeline of work.

#the sed steps refer to formatting steps that get the outputs of the python script ready for use in further python scripts are R scripts used for further analysis.

#Aligning and calling variants for the alternative haplotype. Creating the vcf plus filtered and alterntative chromosome versions
    #Manual
        alt_assembly=
        alt_GCA=
        ref_assembly=
        ref_GCA=

    #Looped sets, create a textfile of your assembiles and GCA numbers to automatically loop through the code.
        big_textfile=path/to/text/file/needed/here
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            alt_GCA=$(echo ${line} | cut -d ";" -f 2 )
            ref_assembly=$(echo ${line} | cut -d ";" -f 3 )
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${alt_GCA} ${ref_assembly} ${ref_GCA}

            #Map the alt to the reference assemblies and make bam file
            output=path/to/output/file
            alternate=path/to/alternative/haplotype/fasta
            ref=path/to/reference/haplotype/fasta
            minimap2 -ax asm10 ${ref} ${alternate} | samtools view -Sb > ${output} 

            output=path/to/output/file/from/minimap2
            samtools_output=path/to/samtools/output/file
            samtools sort ${output} > ${samtools_output}
            samtools index ${samtools_output}

            bcftools mpileup -f ${ref} --min-MQ 60 --annotate FORMAT/DP ${samtools_output} > ${samtools_output}.mpileup
            bgzip ${samtools_output}.mpileup
            bcftools index -f ${samtools_output}.mpileup.gz

            #Turning VCF into bam 
            vcf_output=path/to/bcftools/output/VCF/file/you/are/creating/here
            bcftools call -m --ploidy 1 ${samtools_output}.mpileup.gz | bgzip > ${vcf_output}
            bcftools index ${vcf_output}.gz
            

            #Filtering the VCF to DP=1 
            vcf_depth_filter=path/to/bcftools/filter/output/
            bcftools filter -i 'FORMAT/DP <= 1 & INFO/DP<=1' ${vcf_output} | bgzip > ${vcf_depth_filter}
            bcftools index ${vcf_depth_filter}
    
        done
        
        #Getting the alternative vcf numbered names rather than long name for genome annotation analysis <- REQUIRES MANUAL INPUT TO CHECK RENAME FILE 
        big_textfile=path/to/text/file/needed/here
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            alt_GCA=$(echo ${line} | cut -d ";" -f 2 )
            ref_assembly=$(echo ${line} | cut -d ";" -f 3 )
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${alt_GCA} ${ref_assembly} ${ref_GCA}

            ref=path/to/reference/haplotype/fasta
           
            #get index, get first column
            rename=path/to/rename/text/file/with/long/chrom/names/then/tab/numeric/chrom/names
            samtools faidx ${ref}
                
            #add row number column
            cut -f1 ${ref}.fai | awk '{print $0,NR}' | > ${rename}
            cp ${rename} path/to/create/a/second/version/of/rename/file/that/you/will/alter/to/have/just/the/number/chroms/
            
            Num_text=path/to/create/a/second/version/of/rename/file/that/you/will/alter/to/have/just/the/number/chroms/
            #you have to check that they rename numbers will be correct. If there are some unscaffolded contigs in the middle of the scaffolded ones then the numbers will be wrong!
            #For num text you have to remove all the long chrom names and just leave the numbers of autosomes that you want to analyse in later steps!
            code ${Num_text}
            code ${rename}
        done
    
        #use bcftools to change name in vcf file - Check rename file first!
        big_textfile=path/to/text/file
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            alt_GCA=$(echo ${line} | cut -d ";" -f 2 )
            ref_assembly=$(echo ${line} | cut -d ";" -f 3 )
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${alt_GCA} ${ref_assembly} ${ref_GCA}
            
            vcf_depth_filter=path/to/bcftools/filter/output/
            vcf_depth_filter_chrom_alt=path/to/bcftools/filter/alternative/chromosome/name/output/
            rename=path/to/rename/text/file/with/long/chrom/names/then/tab/numeric/chrom/names
            
            bcftools annotate ${vcf_depth_filter} --rename-chrs ${rename} > ${vcf_depth_filter_chrom_alt}
            bgzip ${vcf_depth_filter_chrom_alt}
            bcftools index -f ${vcf_depth_filter_chrom_alt}.gz
        done

#Getting four fold degenerate sites for whole aligned haplotypes:
        big_textfile=path/to/textfile/
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            alt_GCA=$(echo ${line} | cut -d ";" -f 2 )
            ref_assembly=$(echo ${line} | cut -d ";" -f 3 )
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${alt_GCA} ${ref_assembly} ${ref_GCA}
            asm=ASM10
        
            code_path=path/to/where/scripts/are/stored
            rename=path/to/rename/text/file/with/long/chrom/names/then/tab/numeric/chrom/names
            gff=path/to/genome/annotation/file/from/DTOL
            ref=path/to/reference/haplotype/fasta
            vcf_depth_filter_chrom_alt=path/to/bcftools/filter/alternative/chromosome/name/output/
            four_D_output=path/to/4D/out/put
            four_D_bed_output=path/to/4D/bed/file/output

            python3 ${code_path}codingSiteTypes.py -a ${gff} -f gff3 -r ${ref} -v ${vcf_depth_filter_chrom_alt} --scaffoldLookup ${rename} --ignoreConflicts --useAnnotationScaffoldNames | bgzip > ${four_D_output}
            
            #Turns the output into a bed file format
            gunzip -c  ${four_D_output} | awk 'BEGIN {OFS="\t"}; $5=="4" {print($1"\t"$2-1"\t"$2)}'|  > ${four_D_bed_output}

            #Creates a vcf file of 4D positions only
            output=output/to4D/site/vcf/files
            bedtools intersect -wa -a ${vcf_depth_filter_chrom_alt} -b ${four_D_bed_output} -header > $output 
        done

        #getting 4D snp positions:
        big_textfile=path/to/textfile/
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            set=$(echo ${line} | cut -d ";" -f 2 )
            echo ${alt_assembly} ${set}
            asm=ASM10
            code_path=path/to/where/scripts/are/stored
            vcf_fourD=output/to4D/site/vcf/files
            Assembly_info=path/to/assembly/information/output/
            chrom_text=path/to/textfile/with/list/of/autosome/names/to/be/analysed 

            python3 ${code_path}/Snp_finder.py -i ${vcf_fourD} -i2 ${chrom_text} -o ${Assembly_info} 

            #Gets the output ready for analysis in R
            Assembly_info=path/to/assembly/information/output/
            Assembly_info_sed=path/to/assembly/information/output/sedded
            sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}


#Autosome analysis 
    textfile=path/to/textfile/
    for line in $(cat ${textfile})
        do 
        alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
        echo ${alt_assembly}
   
        code_path=path/to/where/scripts/are/stored
        vcf_depth_filter_chrom_alt=path/to/bcftools/filter/alternative/chromosome/name/output/
        Assembly_info=path/to/assembly/information/output/
        indel_finder=path/to/output/
        indel_binner_100=path/to/output/
        indel_binner_10000=path/to/output/
        chrom_text=path/to/textfile/with/list/of/autosome/names/to/be/analysed 
        every_indel_length=path/to/output/

        python3 ${code_path}/Assembly_information.py -i ${vcf_depth_filter_chrom_alt} -i2 ${chrom_text} -o ${Assembly_info} -o2 ${indel_finder} -o3 ${indel_binner_100} -o4 ${indel_binner_10000} -o5 ${every_indel_length}

        #Getting the outputs ready for use in R
        Assembly_info=path/to/assembly/output/
        Assembly_info_sed=path/to/assembly/output/sedded/output
        sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}

        indel_finder=path/to/output/
        indel_finder_sed=path/to/assembly/output/sedded/output
        sed 's/\[//g' ${indel_finder} | sed 's/\]//g' | > ${indel_finder_sed}


        indel_binner_100=path/to/output/
        indel_binner_100_sed=path/to/assembly/output/sedded/output
        sed 's/{//g' ${indel_binner_100} | sed 's/}//g' | sed 's/:/,/g' | sed 's/'\''//g' | > ${indel_binner_100_sed}
        
        indel_binner_10000=path/to/output/
        indel_binner_10000_sed=path/to/assembly/output/sedded/output
        sed 's/{//g' ${indel_binner_10000} | sed 's/}//g' | sed 's/:/,/g' | sed 's/'\''//g' | > ${indel_binner_10000_sed}

        every_indel_length=path/to/output/
        every_indel_length_sed=path/to/assembly/output/sedded/output
        sed 's/{//g' ${every_indel_length} | sed 's/}//g' | sed 's/:/,/g' | sed 's/'\''//g' | > ${every_indel_length_sed}
    done
    done
#CDS analysis 
    #Pull CDS pos out of gff file:
        big_textfile=path/to/textfile/

        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${ref_GCA}
       
            gff=path/to/genome/annotation/file/from/DTOL
            final_file=path/to/output/
            final_file_short=path/to/output/
            final_file_sorted=path/to/output/
            final_file_sorted_merged=path/to/output/
            final_file_sorted_short=path/to/output/
            final_file_sorted_merged_short=path/to/output/
            #pull out all the cds rows
                awk '{if($3=="CDS") print $0}' ${gff} > all_exons.bed
            #remove and append column 1 and 9
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
        #Making the CDS included file for the vcf with the number chrom 
            textfile=path/to/textfile/
            for line in $(cat ${textfile})
                do 
                alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                set=$(echo ${line} | cut -d ";" -f 2 )
                echo ${alt_assembly} ${set}

        
                CDS_file=path/to/file/made/in/previous/step/with/gff/info
                vcf_depth_filter_chrom_alt=path/to/bcftools/filter/alternative/chromosome/name/output/
                CDS_include=path/to/cds/included/vcf/file/

                bcftools view ${vcf_depth_filter_chrom_alt} --targets-file ${exon_file} -o ${output_include}
                bgzip -f ${output_include}
                bcftools index ${output_include}.gz
                
                #Getting snps for whole cds regions - need for indel diversity denominator 
                code_path=path/to/where/scripts/are/stored
                CDS_include=path/to/cds/included/vcf/file/
                Assembly_info=path/to/assembly/information/output/
                chrom_text=path/to/textfile/with/list/of/autosome/names/to/be/analysed 
            
                python3 ${code_path}/Snp_finder.py -i ${output_include}.gz -i2 ${chrom_text} -o ${Assembly_info} 

                Assembly_info=path/to/assembly/information/output/
                Assembly_info_sed=path/to/output/sed/output
                  sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}


                #getting four d positions in CDS and then snps at the four d positions 
                four_D_bed_output=output/to4D/site/bed/file
                output=path/to/output/
                bedtools intersect -wa -a ${output_include}.gz -b ${four_D_bed_output} -header > $output 
                

                code_path=path/to/where/scripts/are/stored
                vcf_fourD=output/to4D/site/vcf/files
                Assembly_info=path/to/assembly/information/output/
                chrom_text=path/to/textfile/with/list/of/autosome/names/to/be/analysed 
            
                
                python3 ${code_path}/Snp_finder.py -i ${vcf_fourD} -i2 ${chrom_text} -o ${Assembly_info} 


                Assembly_info=path/to/assembly/information/output/
                Assembly_info_sed=path/to/assembly/output/sedded/output
                sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}
    #pulling assembly information for CDS region
        textfile=path/to/textfile/
        for line in $(cat ${textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            echo ${alt_assembly}
        
            CDS_file=file/path/to/cds/positions
            vcf_depth_filter_chrom_alt=path/to/bcftools/filter/alternative/chromosome/name/output/
            CDS_include=path/to/cds/included/vcf/file/

            code_path=path/to/where/scripts/are/stored
            CDS_include=path/to/cds/included/vcf/file/
            Assembly_info=path/to/output/
            indel_finder=path/to/output/
            indel_binner_100=path/to/output/
            indel_binner_10000=path/to/output/
            chrom_text=path/to/textfile/with/list/of/autosome/names/to/be/analysed 

            python3 ${code_path}Assembly_information.py -i ${CDS_include} -i2 ${chrom_text} -o ${Assembly_info} -o2 ${indel_finder} -o3 ${indel_binner_100} -o4 ${Indel_binner_10000}


        #Sorting of the outputs for R 
            #sorting indel pos output - can then be loaded into R as a table before R then calculates the length of them! 
                sed_indel_pos_input=path/to/output/
                indel_pos_sed_output=path/to/assembly/output/sedded/output
                sed 's/\[//g' ${sed_indel_pos_input} |  sed 's/\]//g' |sed 's/\,//g' | sed "s/'//g" | > ${indel_pos_sed_output}
            
            #Sorting indel finder output - just removes brackets and stuff for loading into R
                input_sed_indelfinder=path/to/output/
                IndelFinder_sed_output=path/to/assembly/output/sedded/output
                sed 's/\[//g' ${input_sed_indelfinder} | sed 's/\]//g' | head > ${IndelFinder_sed_output}
                
            #Sorting <100 indel binning output - removes brackets, gets column 1 to all the same decimal places so they can be sorted correctly and removes empty columns - also removes the > and < columns but these can still be seen in the orginal output
                indel_binner_output=path/to/output/
                indel_binner_output_sed=path/to/assembly/output/sedded/output
                indel_binner_output_sed_r=path/to/assembly/output/sedded/output/R

                sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > save.txt
                awk '{print $1}' save.txt |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                awk '{print $2}' save.txt | > col_2.txt
                paste col_1.txt col_2.txt > unsorted.txt
                awk '$1!=0000' unsorted.txt > sorted.txt
                sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r} 

            #sorting indel binner up to 10,000 - more complex than <100 output as column 1 needs more work but does the same thing pretty much - also removes the > and < columns but these can still be seen in the orginal output 
                indel_binner_output_10000=path/to/output/
                indel_binner_output_10000_sed=path/to/assembly/output/sedded/output
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
            big_textfile=path/to/textfile/
            for line in $(cat ${big_textfile})
                do 
                alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                echo ${alt_assembly}

                #Getting bed files rather than text, via renaming and then altering of file to get proper format
                    indel_pos_sed_output_text=path/to/output/
                    rename_text_sed2=path/to/output/with/new/bed/ending
                    cp ${indel_pos_sed_output_text} ${rename_text_sed2}
                
                    CDS_file_short=path/to/CDS/postions/file
                    CDS_file_short_sed=path/to/output/sed
                    indel_pos_sed_output=path/to/output/
                    indel_pos_sed_output_sed=path/to/assembly/output/sedded/output
                    output=path/to/output/
                
                #Bedtools hates spaces so have to makesure input is tabbed and also that the indel pos jave 0 based column 2 and 1 based column 3:
                    sed 's/ /\t/g' ${exon_file_short} > ${exon_file_short_sed}
                    awk '{print $2-1}' $rename_text_sed2 > column_2.txt
                    awk '{print $1}' $rename_text_sed2 > column_1.txt
                    awk '{print $3}' $rename_text_sed2 > column_3.txt
                    paste column_1.txt column_2.txt column_3.txt > ${indel_pos_sed_output}
                    sed 's/ /\t/g' ${indel_pos_sed_output} | sed 's/,//g' |sed 's/'\''//g' | > ${indel_pos_sed_output_sed}

                #Actual intersect command 
                    bedtools intersect -wao -a ${indel_pos_sed_output_sed} -b ${exon_file_short_sed} > ${output}

        #R-STEPS - in R you sort the bedtools intersect output so that it is split into in cds and outwith cds region indels. We ignore indels that aren't fully in either. Indel lengths are also got in this step by substracting the start indel positions from end positions in R. 
            #this gets the assemblies ready to go through the indel dictionary sorter from r
            #Indels in CDS:
                text_file=path/to/textfile
                for line in $(cat ${text_file})
                    do 
                    alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                    echo ${alt_assembly} 

                    indels_in_CDS_length=path/to/output/from/R
                    output=path/to/output/
                    
                    awk 'NR!=1{print $1}' ${indels_in_CDS_length} > ${output}

                #Indel dictionary sorter from R
                    text_file=path/to/textfile
                    for line in $(cat ${text_file})
                        do 
                        alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                        echo ${alt_assembly} 
                        code_path=path/to/where/scripts/are/stored
                        input=path/to/output/from/r/of/indel/lengths
                        dict_100=path/to/output/
                        greater_than_100_dict=path/to/output/

                        python3 ${code_path}/Indel_sorter.py -i ${input} -o ${dict_100} -o2 ${greater_than_100_dict}

                        indel_binner_output=/path/to/output/
                        indel_binner_output_sed=path/to/assembly/output/sedded/output
                        indel_binner_output_sed_r=path/to/assembly/output/sedded/output/R

                        sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > save.txt
                        awk '{print $1}' save.txt |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                        awk '{print $2}' save.txt | > col_2.txt
                        paste col_1.txt col_2.txt > unsorted.txt
                        awk '$1!=0000' unsorted.txt > sorted.txt
                        sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                        sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r} 

                #Indels outside CDS 
                text_file=path/to/textfile
                    for line in $(cat ${text_file})
                        do 
                        alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                        echo ${alt_assembly} 
                
                input_file=path/to/output/of/indels/outwith/cds/region/from/r
                input_file_2=path/to/output/sedded
                  sed 's/" //g' $input_file | sed 's/"//g' | sed 's/ //g' | sed 's/\t/,/g' | > $input_file_2


                    text_file=path/to/textfile
                    for line in $(cat ${text_file})
                        do 
                        alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                        echo ${alt_assembly} 
                        code_path=path/to/where/scripts/are/stored
                        input_file_2=/path/to/output/sedded/from/above
                        dict_100=/path/to/output/
                        greater_than_100_dict=path/to/output
                     

                        python3 ${code_path}/Indels_sorter.py  -i ${input_file_2} -o $ ${dict_100} -o2 ${greater_than_100_dict} 

                        indel_binner_output=path/to/output/
                        indel_binner_output_sed=path/to/assembly/output/sedded/output
                        indel_binner_output_sed_r=path/to/assembly/output/sedded/output/R
                        sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > save.txt
                        awk '{print $1}' save.txt |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                        awk '{print $2}' save.txt | > col_2.txt
                        paste col_1.txt col_2.txt > unsorted.txt
                        awk '$1!=0000' unsorted.txt > sorted.txt
                        sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                        sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r}
                        

#Repeats analysis assembly info
    #Repeat_Regions, runs through the reference fasta and pulls out the positions of repeat regions
        big_textfile=path/to/textfile/
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            ref_assembly=$(echo ${line} | cut -d ";" -f 3)
            ref_GCA=$(echo ${line} | cut -d ";" -f 4 )
            echo ${alt_assembly} ${ref_assembly} ${ref_GCA}

            code_path=path/to/where/scripts/are/stored
            chrom_names=path/to/textfile/with/autosome/names
            ref=path/to/reference/haplotype/fasta
            Repeat_positions=path/to/output/
            Repeat_positions_sed=path/to/output/sedded
            #repeat regions come out with correct bed file formatiing (col2=0 col3=1 based)

            python3 ${code_path}/Repeat_finder.py -i ${ref} -i2 ${chrom_names_R} -o ${Repeat_positions} 
            sed 's/\[//g' ${Repeat_positions} | sed 's/\]//g' | > ${Repeat_positions_sed}
  
    #Intersecting the repeat regions with indels not in cds
        big_textfile=path/to/textfile/
        for line in $(cat ${big_textfile})
            do 
            alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
            echo ${alt_assembly} 
            
            indels_not_in_CDS=path/to/output/
            indels_not_in_CDS_sed=path/to/output/sedded
            awk 'NR!=1{print $0}' ${indels_not_in_CDS} | sed 's/ //g' | sed 's/"//g' | > ${indels_not_in_CDS_sed}
          
            Repeat_positions_sed=path/to/output/sedded
            repeat_pos_ready=path/to/output/sedded/again
            seddy_output=path/to/assembly/output/sedded/output
            repeat_pos_final=path/to/output/final/for/this/section
            rename=path/to/rename/text/file/with/long/chrom/names/then/tab/numeric/chrom/names
            
            sed 's/'\''//g' $Repeat_positions_sed | sed 's/,//g' | > ${repeat_pos_ready}

            join ${repeat_pos_ready} ${rename} > ${seddy_output}
            awk '{print $4"\t"$2"\t"$3}' $seddy_output | > ${repeat_pos_final}
           
    
            output=path/to/output/
            bedtools intersect -wao -a ${indels_not_in_CDS_sed} -b ${repeat_pos_final} > ${output}

        #Python analysis - in repeats - go to R and split up the indels into in repeats and outwith repeats. Indel lengths are also got in this step by substracting the start indel positions from end positions in R. 
        text_file=path/to/textfile
                for line in $(cat ${text_file})
                    do 
                    alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                    echo ${alt_assembly} 

                    input=path/to/output/ofindel/pos/in/repeats
                    output=path/to/output/of/indel/pos/in/repeats/sedded

                    awk 'NR!=1{print $0}' ${input} | sed 's/ //g' | sed 's/""/"/g' | sed 's/"/,/g' | sed 's/,$//g' | > ${output}
                    
                text_file=path/to/textfile
                for line in $(cat ${text_file})
                    do 
                    alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                    echo ${alt_assembly} 

                    code_path=path/to/where/scripts/are/stored
                    input_file=path/to/output/of/repeat/positions
                    dict_100=/path/to/output/
                    greater_than_100_dict=path/to/output/
              

                    python3 ${code_path}/Indels_sorter.py -i ${input_file_2} -o1 ${dict_100} -o1 ${greater_than_100_dict} 

                    indel_binner_output=path/to/output/
                    indel_binner_output_sed1=path/to/assembly/output/sedded/output1
                    indel_binner_output_sed=path/to/assembly/output/sedded/output
                    indel_binner_output_sed_r=path/to/assembly/output/sedded/output/R
                    sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > ${indel_binner_output_sed1}
                    awk '{print $1}' ${indel_binner_output_sed1} |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                    awk '{print $2}' ${indel_binner_output_sed1} | > col_2.txt
                    paste col_1.txt col_2.txt > unsorted.txt
                    awk '$1!=0000' unsorted.txt > sorted.txt
                    sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                    sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r}
    

        #Python analysis - out repeats

                text_file=path/to/textfile
                for line in $(cat ${text_file})
                    do 
                    alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                    echo ${alt_assembly} 

                    input=path/to/output/of/indels/not/in/repeats
                    output=path/to/output/indels/not/in/repeats/sedded

                    awk 'NR!=1{print $0}' ${input} | sed 's/ //g' | sed 's/""/"/g' | sed 's/"/,/g' | sed 's/,$//g' | > ${output}
                    
                    
                    text_file=path/to/textfile
                    for line in $(cat ${text_file})
                        do 
                        alt_assembly=$(echo ${line} | cut -d ";" -f 1 )
                        echo ${alt_assembly} 
                    
                        code_path=path/to/where/scripts/are/stored
                        input_file=path/to/output/of/indels/not/in/repeats/sedded
                        indel_total_lengths=path/to/output/
                        dict_100=path/to/output/
                        greater_than_100_dict=path/to/output/
                        Number_of_indels=path/to/output/

                        python3 ${code_path}/Indels_sorter.py -i ${input_file} -o ${dict_100} -o2 ${greater_than_100_dict} 

                        indel_binner_output=path/to/output/
                        indel_binner_output_sed1=path/to/assembly/output/sedded/output1
                        indel_binner_output_sed=path/to/assembly/output/sedded/output
                        indel_binner_output_sed_r=path/to/assembly/output/sedded/output/r
                        sed 's/{//g' ${indel_binner_output} | sed 's/}//g' |  sed 's/'\''//g' |  sed 's/:/'\ '/g' | > ${indel_binner_output_sed1}
                        awk '{print $1}' ${indel_binner_output_sed1} |  awk '{printf ("%03d\n", $1)}' | > col_1.txt
                        awk '{print $2}' ${indel_binner_output_sed1} | > col_2.txt
                        paste col_1.txt col_2.txt > unsorted.txt
                        awk '$1!=0000' unsorted.txt > sorted.txt
                        sort -k 1 sorted.txt > ${indel_binner_output_sed} 
                        sed 's/ /\t/g' ${indel_binner_output_sed} | sed 's/\t/,/g' | sed 's/$/,/g' | > ${indel_binner_output_sed_r}
                    

                #making vcfs from the repeat regions
                Repeat_positions_vcf=path/to/output/
                vcf_depth_filter_chrom_alt=path/to/output/
                repeats_output=path/to/output/

                bcftools view ${vcf_depth_filter_chrom_alt} --targets-file ${Repeat_positions_vcf} -o ${repeats_output}
                bgzip -f ${repeats_output}
                bcftools index ${repeats_output}.gz


    #Getting snps for whole repeat regions - need for indel diversity denominator, already have indels from intersect!
    code_path=path/to/where/scripts/are/stored
    repeats_output=path/to/output/
    Assembly_info=path/to/assembly/information/output/
    chrom_text=path/to/textfile/with/list/of/autosome/names/to/be/analysed 

    python3 ${code_path}/Snp_finder.py -i ${repeats_output}.gz -i2 ${chrom_text} -o ${Assembly_info} 

    Assembly_info=path/to/assembly/information/output/
    Assembly_info_sed=path/to/assembly/output/sedded/output
    sed 's/\[//g' ${Assembly_info} | sed 's/\]//g' | > ${Assembly_info_sed}
