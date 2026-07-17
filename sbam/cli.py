import argparse
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(
        prog="sbam",
        description="SBAM: Structural and Base-Level analysis of ONT microbial assemblies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required inputs
    parser.add_argument("-a", "--assembly", required=True, 
                        help="Path to the consensus assembly (FASTA)")
    parser.add_argument("-r", "--reads", required=True, 
                        help="Path to the raw long reads (FASTQ)")
    parser.add_argument("-o", "--outdir", required=True, 
                        help="Directory to save the BAM files and HTML report")
    
    # Optional parameters
    parser.add_argument("-t", "--threads", type=int, default=4, 
                        help="Number of CPU threads for minimap2")
    parser.add_argument("--buffer-size", type=int, default=50000, 
                        help="Size of the cyclic buffer edge (bp) for mapping")
    parser.add_argument("--read-type", choices=["map-ont", "map-pb"], default="map-ont", 
                        help="Sequencing technology used (default: map-ont)")
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    # Setup output directory
    os.makedirs(args.outdir, exist_ok=True)
    print(f"[INIT] Starting SBAM analysis...")
    print(f"[INIT] Assembly: {args.assembly}")
    print(f"[INIT] Output directory: {args.outdir}\n")
  
# Structural validation step
    print("Structural Validation")

    from sbam.struct.struct import GenomeArchitecture 
    from sbam.struct.mapping import ReadAligner
    from sbam.struct.eval import JunctionEvaluator
    
    # Initialize the architecture class
    genome = GenomeArchitecture(args.assembly, args.outdir)
    print(f" > Found {len(genome.records)} contig(s).")
    
    # Create the cyclic buffer
    cyclic_fasta, orig_lengths = genome.create_cyclic_buffer(buffer_size=args.buffer_size)
    print(f" > Cyclic buffer applied. Saved to: {cyclic_fasta}")
        
    # Align raw reads to the cyclic buffer
    aligner = ReadAligner(args.reads, cyclic_fasta, args.outdir, args.threads)
    bam_path = aligner.run_mapping(read_type=args.read_type)
    
    # Evaluate Junctions
    evaluator = JunctionEvaluator(bam_path, orig_lengths)
    junction_metrics = evaluator.evaluate_junctions()
    
    print("\n Junction Structural Scores")
    for contig, metrics in junction_metrics.items():
        length_kb = metrics['length'] / 1000
        depth = metrics['avg_depth']
        
        # Format the output
        print(f"   - {contig} ({length_kb:.1f} kb | Depth: {depth}x): "
              f"Score {metrics['spanning_score']:.2f} [{metrics['status']}] "
              f"({metrics['spanning_reads']} spanning / {metrics['broken_reads']} broken)")

    from sbam.struct.ori import BioVal
    
    physics = BioVal(args.assembly)
    physics_metrics = physics.calculate_symmetry()
    
    print("\n Biological scores")
    for contig, metrics in physics_metrics.items():
        print(f"   - {contig}: Symmetry {metrics['symmetry']}° [{metrics['viability']}] "
              f"(oriC ≈ {metrics['oric_pos']} bp, ter ≈ {metrics['ter_pos']} bp)")    
        
    # Biological Validation
    print("\n Base-Level Validation")

    from sbam.base.base import BaseAccuracy
    
    # Calculate base-level accuracy and identify systematic error motifs
    fidelity = BaseAccuracy(bam_path=bam_path, fasta_path=args.assembly, threads=args.threads)
    
    # analyze_motifs() automatically builds the profile, masks structural noise, 
    # and returns a list of statistically significant discordant motifs.
    motif_results = fidelity.analyze_motifs()
    
    print("\n > --- Top Systematic Error Motifs (by k-mer size) ---")
    
    if not motif_results:
        print("   No significant systematic error motifs detected. Basecalling is highly concordant.")
    else:
        # Group the results by k-mer size for clean printing
        for k in [4, 5, 6]:
            print(f"\n   [{k}-mers]")
            # Filter results for this specific k-mer length
            k_results = [res for res in motif_results if res['kmer_size'] == k]
            
            if not k_results:
                print("     No significant motifs detected.")
            else:
                for res in k_results[:3]: # Print the top 3 motifs for each k-mer size
                    print(f"     {res['motif']} | Discordance Rate: {res['error_rate']:.1%} "
                          f"({res['discordant_count']}/{res['occurrences']} hits) | "
                          f"Fold Increase: {res['fold_increase']:.1f}x")

    # Reporting
    print("\n[PHASE 3] Generating Dashboard")
    from sbam.report.report import DashboardBuilder
    
    builder = DashboardBuilder(args.outdir, args.assembly)
    masked_pct = fidelity.masked_pct if hasattr(fidelity, 'masked_pct') else 0
    masked_bases = fidelity.masked_bases if hasattr(fidelity, 'masked_bases') else 0
    
    report_file = builder.generate_report(
        junction_metrics=junction_metrics,
        physics_metrics=physics_metrics,
        motif_results=motif_results,
        masked_bases=masked_bases, 
        masked_pct=masked_pct
    )
    
    print(f"\n[DONE] SBAM execution finished successfully.")
    print(f"       View your report at: {report_file}")
    return 0