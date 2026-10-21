import sys, os
import main, tools, cui, seq_io

########################################################################
class MetaSorter:
    def __init__(self, version="", date_of_creation=""):
        self.version = version
        self.date_of_creation = date_of_creation
    
        # Input/output
        self.oIO = seq_io.IO()
        

    def parse_gff_file(self, gff_path, contigs):
        # Parse GFF file
        if not gff_path or not os.path.exists(gff_path):
            raise ValueError(f"Input GFF file {gff_path} does not exists!")
        tools.msg(f"Reading {gff_path}")
        dGFF = self.oIO.readGFF(path=gff_path, contigs=contigs, mode='dictionary')
        dGFF['input file'] = gff_path
        return [dGFF]
        
    def select_bin_gff(self, master_gff, contigs):
        records = [rec for rec in master_gff if rec.startswith(tuple(contigs))]
        return records
    
    def execute(self, args):       
        # Set runtime options
        options = cui.get_options()
        options['-u'] = args.input_folder               # Input folder
        options['-o'] = args.output_folder              # Output folder
        options['-d'] = args.project_folder             # Project folder
        options['-d'] = args.project_folder             # Project folder
        options['-c'] = args.mqs                        # Minimum methylation/modification quality score
        options['-ogf'] = args.output_graph_format      # Graphical format
        options['-dpi'] = args.dpi                      # DPI for rastr format
        options['-wl'] = args.sliding_window_length     # Sliding window length
        options['-ws'] = args.sliding_window_step       # Sliding window step
        options['-mm'] = args.circular_graph.upper()    # Request circular graph
        options['-sp'] = "Y"                            # Request statistic panel
        options['-s'] = "sites"                         # Request site statistics instead of motif statistics
        options['-f'] = "M"                             # Request statistics of methylated sites instead of unmethylated 'U'
        options['-r'] = 0                               # No limits for number of validated sites (server application)
        options['-x'] = os.path.join("lib","bin") if args.bin_path == "" else args.bin_path         # Path to directory with executable blast files 
        options['-tmp'] = os.path.join("lib","bin","tmp") if args.tmp_path == "" else args.tmp_path # Path to temporary blast directory s 
        
        separator = args.file_name_separator
        
        # List of motifs
        motifs = [s.strip() for s in args.motifs.split(";")]
        
        if len(motifs) == 0:
            sys.exit("No motifs were specified!")
            
        # -------------------------------------------------------
        # Validate grouping parameters
        # -------------------------------------------------------
        
        if not 0 < args.grouping_p_cutoff < 1:
            sys.exit("ERROR: --grouping_p_cutoff must be between 0 and 1 (exclusive)!")
        
        if args.min_motifs < 1:
            sys.exit("ERROR: --min_motifs must be a positive integer (>= 1)!")
        
        if args.min_coverage < 1:
            sys.exit("ERROR: --min_coverage must be a positive integer (>= 1)!")
        
        if not 0 < args.alpha < 1:
            sys.exit("ERROR: --alpha must be between 0 and 1 (exclusive)!")
        
        if not 0 <= args.max_group_difference <= 1:
            sys.exit("ERROR: --max_group_difference must be between 0 and 1 (inclusive)!")
        
        if not 0 <= args.consensus_threshold <= 1:
            sys.exit("ERROR: --consensus_threshold must be between 0 and 1 (inclusive)!")
        
        if not 0 <= args.min_motif_weight <= 1:
            sys.exit("ERROR: --min_motif_weight must be between 0 and 1 (inclusive)!")
        
        if not 0 <= args.max_motif_weight <= 1:
            sys.exit("ERROR: --max_motif_weight must be between 0 and 1 (inclusive)!")
        
        if args.min_motif_weight > args.max_motif_weight:
            sys.exit("ERROR: --min_motif_weight cannot exceed --max_motif_weight!")
                
        fasta_extensions = (".FA", "FASTA", ".FNA", ".FST")
        genbank_extensions = (".GB", ".GBK", ".GBF", ".GBFF")
        # Current version works only with FASTA files
        genbank_extensions = ()
        
        input_path = os.path.join(args.input_folder, args.project_folder)
        master_gff = None
        if args.gff_file:
            master_gff_file_path = os.path.join(input_path, args.gff_file)
            if os.path.exists(master_gff_file_path):
                master_gff = self.oIO.open_text_file(path = master_gff_file_path, flg_inlist = True)
        
        # Loop through pairs of GFF and sequence files
        for seq_fname in [fn for fn in os.listdir(os.path.abspath(os.path.join(args.input_folder, args.project_folder))) 
                if fn.upper().endswith(fasta_extensions + genbank_extensions)]:
            basename = seq_fname.split(separator)[0] if separator else seq_fname[:seq_fname.rfind(".")]
            outpath = os.path.join(args.output_folder, args.project_folder, basename)
            os.makedirs(outpath, exist_ok=True)
            
            seq_path = os.path.join(input_path, seq_fname)
            # Parse sequence file
            genome = self.oIO.get_SeqDataset(seq_path, fname_separator = separator)
            contigs = genome.get_contigs()
            contig_locations = genome.get_contig_locations()
            contig_titles = genome.get_contig_titles(flg_with_locations=True)
                
            if not genome or not contigs:
                raise ValueError(f"File {seq_fname} has wrong format or corrupted!")
                    
           
            # Parse GFF file    
            if master_gff is None:
                candidate_gff_files = [fn for fn in os.listdir(os.path.abspath(os.path.join(args.input_folder, args.project_folder))) 
                    if fn.upper().startswith(basename.upper()) and fn.upper().endswith(".GFF")]
                
                # Skip if pair sequence file has not been found
                if len(candidate_gff_files) == 0:
                    tools.alert(f"GFF file for sequence file {seq_fname} in folder {os.path.join(args.input_folder, args.project_folder)} has not been found!")
                    continue
                    
                   
                gff_fname = candidate_gff_files[0]
                gff_path = os.path.join(input_path, gff_fname)
                
            else:
                gff_fname = basename + ".gff"
                gff_path = os.path.join(outpath, gff_fname)
                gff_records = self.select_bin_gff(master_gff, genome.get_contig_titles())
                if len(gff_records) == 0:
                    tools.alert(f"GFF records for sequence file {seq_fname} have not been found!")
                    continue
                    
                self.oIO.save("\n".join(gff_records), gff_path)
            
            
            GFF_list = self.parse_gff_file(gff_path, contigs)
            
            
            # Groups of contigs
            group_collection = []
            # filtered contigs
            filtered_contigs = []
            try:
                filter_contigs = float(args.filter_chimeric_contigs)
                if filter_contigs > 0.05 or filter_contigs <= 0:
                    filter_contigs = False
            except:
                filter_contigs = True if args.filter_chimeric_contigs.upper()[0] == "Y" else False
            
            # Loop through motifs
            
            for motif in motifs:
                options['-w'] = motif
                tools.msg(f"Searching for motif {motif} in files {gff_fname} and {seq_fname}...")
                # Run program

                oMain = main.Interface(options, self.version, self.date_of_creation)
                outputs = oMain.execute(genome, contigs, GFF_list)
                
                if not outputs:
                    continue
    
                for output in outputs:
                    # groups = {'motif':motif, 'positive':[], 'negative':[], 'median':[], 'filtered_contigs':[], 'grouped':False}
                    groups = oMain.group_contigs(output=output, 
                        motif=motif, 
                        p_cutoff=args.grouping_p_cutoff,                 # Cutoff for biased distribution within contigs
                        min_motifs=args.min_motifs,                      # Minimal number of modified sites per contig
                        min_coverage=args.min_coverage,                  # Minimal average coverage at modified sites
                        alpha=args.alpha,                                # alpha cutoff for modified motif distribution across motifs 
                        max_group_difference=args.max_group_difference,  # Minimal difference between groups in fractions of modified motifs
                        filter_contigs=filter_contigs,
                        )
                    filtered_contigs += groups['filtered_contigs']
        
                    if groups['grouped'] is False:
                        continue
                    
                    # Append grouped contigs
                    group_collection.append(groups)
                    
                    # Remove graphs from outputs if not needed
                    if args.save_graphs.upper() == "N" and 'svg' in output:
                        output['svg'] = None
                        
                    # Save output
                    oMain.save_output(output, path=outpath)
                '''
                print("\nmetasorter:192")
                print(groups)
                print()
                '''
                # Clean TMP folder
                oMain.clean_tmp_folder()
                tools.msg("\tDone.")

            '''
            print("metasorter:201", len(group_collection))
            for group in group_collection:
                print()
                print(group)
            print()
            '''

            grouped_contigs = tools.group_contigs(contig_locations, 
                group_collection, 
                unresolved_contigs = args.unresolved_contigs,
                consensus_threshold = args.consensus_threshold, # Minimum weighted agreement (0–1) for grouping contigs.
                                                                # Higher values produce more groups (stricter clustering).
                                                                # Lower values tolerate disagreements between motifs.
                min_motif_weight = args.min_motif_weight,       # Minimum weight assigned to an informative motif.
                                                                # Prevents motifs with weak statistical evidence
                                                                # from having negligible influence on consensus.
                max_motif_weight = args.max_motif_weight,       # Maximum weight assigned to an informative motif.
                                                                # Limits the influence of motifs with very small p-values.
                                                                # Prevents one motif from dominating the consensus.
                max_group_difference=args.max_group_difference, # Minimal difference between groups in fractions of modified motifs 
                                                                # Used to compare low-methylated contigs with unmethylated contigs
                                                                # joined into 'uncertain' group.
                uncertain_weight=args.uncertain_weight,         # Relative weight of evidence involving uncertain contigs.
                                                                # A value of 0.1 means that comparisons involving
                                                                # uncertain contigs contribute 10% of the normalmotif weight.
                                                                # A value of 0 disables uncertain-contig evidence.
                )
            if len(grouped_contigs) == 0:
                tools.msg(f"Failed with {basename}!")
                continue
            
            seq = genome.get_sequence()
            for i in range(len(grouped_contigs)):
                contig_group = grouped_contigs[i]
                outfile_name = f"{basename}.{i + 1}.fasta" if len(grouped_contigs) > 1 else f"{basename}.fasta"
                group = []
                for j in range(len(contig_group)):
                    contig_location = contig_group[j]
                    contig = genome.get_contig_by_location(contig_location)
                    if contig is None:
                        sys.exit(f"ERROR: contig [{contig_location}] was not found!")
                    contig_title = str(contig.qualifiers['name'][0] if 'name' in contig.qualifiers else
                        (contig.qualifiers['id'][0] if 'id' in contig.qualifiers else f"contig.{i + 1}.{j + 1}"))
                        
                    # Filter contigs
                    if contig_title in filtered_contigs:
                        continue
                    group.append(f">{contig_title}\n{seq[contig.location.start:contig.location.end]}")
                # Save non-empty group to FASTA file
                if len(group):
                    self.oIO.save(data = "\n".join(group),
                        path = os.path.join(outpath, outfile_name),
                        data_format = "text"
                    )
                    
            # Save filtered contigs
            if filtered_contigs:
                self.oIO.save(data = "\n".join(filtered_contigs),
                    path = os.path.join(outpath, "filtered_contigs.txt"),
                    data_format = "text"
                )
