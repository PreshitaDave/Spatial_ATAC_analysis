
* To organize everything, I create one object for every tissue (488B and 489) and every depth of sequencing (deepseq and lowseq). 
* Create **barcodes (before and after edge effect filtering), bam files (with all barcodes), and fragment files from the bam files** for each object. Save this in the correct folders and the process to do this in one script. 

We know that "edge effects" exist where spots on the left and the right side have higher nFrags and that can maybe skew the analysis so we get rid of them. But we want to keep track of the barcodes before and after edge effect filtering so we can compare them and see how many barcodes we lose from that step.

