#macs (peak calling) example
macs3 callpeak -t minimap2/results/zbtb24_chip_rep1.bam -c minimap2/results/zbtb24_input_rep1.bam -n zbtb24_rep1 --outdir macs3_output -f BAMPE

#before idr the peaks should be sorted by -log10(p-value) 
sort -k8,8nr zbtb24_rep1_peaks.narrowPeak > sorted/zbtb24_rep1_peaks.narrowPeak 

#a note: before idr, macs should be done less stringent for example adding "-p 1e-3" to macs callpeak 

#idr example
idr --samples zbtb24_rep1_peaks_sorted.narrowPeak zbtb24_rep2_peaks_sorted.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file-type bed \ 
--output-file idr/zbtb24_idr.bed \
--plot \
--log-output-file idr/zbtb24_idr.log

#to only get the merged part from idr peak files:
awk '{if($5 >= 540) print $0}' zbtb24_idr > zbtb24_idr_merged




