import os
import json
from datetime import datetime
from jinja2 import Template

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>SBAM Assembly Report</title>
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: #f4f4f9; color: #333; margin: 0; padding: 20px; line-height: 1.5; }
        .container { max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
        h1 { border-bottom: 2px solid #2c3e50; padding-bottom: 10px; color: #2c3e50; margin-bottom: 5px; }
        h2 { color: #34495e; margin-top: 40px; border-bottom: 1px solid #eee; padding-bottom: 5px; }
        .subtitle { color: #7f8c8d; font-size: 14px; margin-bottom: 30px; }
        
        .executive-summary { background: {{ result_color }}22; border: 2px solid {{ result_color }}; padding: 20px; border-radius: 8px; margin-bottom: 30px; }
        .executive-summary h2 { margin: 0 0 10px 0; color: {{ result_color }}; border: none; padding: 0; }
        .executive-summary p { margin: 0; font-size: 16px; font-weight: 500; color: #2c3e50; }
        
        .grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin-bottom: 30px; }
        .card { padding: 20px; border-radius: 6px; background: #ecf0f1; border-left: 5px solid #3498db; }
        .card h3 { margin: 0 0 5px 0; font-size: 13px; text-transform: uppercase; color: #7f8c8d; }
        .card .value { font-size: 20px; font-weight: bold; color: #2c3e50; }
        .card .desc { font-size: 12px; color: #555; margin-top: 8px; }
        
        table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #f8f9fa; color: #2c3e50; }
        
        .badge-pass, .badge-acceptable { background: #2ecc71; color: white; padding: 4px 8px; border-radius: 12px; font-size: 12px; font-weight: bold; }
        .badge-fail, .badge-warning { background: #e74c3c; color: white; padding: 4px 8px; border-radius: 12px; font-size: 12px; font-weight: bold; }
        .badge-warn, .badge-atypical, .badge-no_data { background: #f1c40f; color: black; padding: 4px 8px; border-radius: 12px; font-size: 12px; font-weight: bold; }
        
        .plot-container { background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 10px; margin-top: 20px; }

        /* Dynamic Explanation Box Styles */
        .explanation-box-success { background: #eafaf1; border: 1px solid #2ecc71; border-radius: 6px; padding: 15px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.02); }
        .explanation-box-success h4 { margin: 0 0 10px 0; color: #27ae60; font-size: 14px; text-transform: uppercase; }
        .explanation-box-success p { margin: 0; font-size: 13px; color: #2c3e50; }

        .explanation-box-warning { background: #fef5e7; border: 1px solid #e67e22; border-radius: 6px; padding: 15px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.02); }
        .explanation-box-warning h4 { margin: 0 0 10px 0; color: #d35400; font-size: 14px; text-transform: uppercase; }
        .explanation-box-warning p { margin: 0; font-size: 13px; color: #2c3e50; }

        /* Tooltip and Link Styles */
        .tooltip { position: relative; display: inline-block; cursor: help; color: #7f8c8d; margin-left: 8px; font-size: 16px; vertical-align: middle; }
        .tooltip .tooltiptext { visibility: hidden; width: 350px; background-color: #34495e; color: #fff; text-align: left; border-radius: 6px; padding: 12px; position: absolute; z-index: 100; bottom: 125%; left: 50%; margin-left: -175px; opacity: 0; transition: opacity 0.2s; font-size: 13px; font-weight: normal; box-shadow: 0px 4px 6px rgba(0,0,0,0.3); line-height: 1.4; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; text-transform: none; }
        .tooltip .tooltiptext::after { content: ""; position: absolute; top: 100%; left: 50%; margin-left: -6px; border-width: 6px; border-style: solid; border-color: #34495e transparent transparent transparent; }
        .tooltip:hover .tooltiptext { visibility: visible; opacity: 1; }
        .docs-link { font-size: 12px; margin-left: 15px; color: #3498db; text-decoration: none; font-weight: normal; vertical-align: middle; }
        .docs-link:hover { text-decoration: underline; }
    </style>
</head>
<body>
    <div class="container">
        <h1>SBAM Assembly Report</h1>
        <div class="subtitle">Assembly: <strong>{{ assembly_name }}</strong> | Generated: {{ date }}</div>
        
        <!-- Executive Summary -->
        <div class="executive-summary">
            <h2>Result: {{ result }}</h2>
            <p>{{ result_msg }}</p>
            
            <h3 style="margin: 25px 0 10px 0; font-size: 15px; color: #2c3e50;">Assembly Overview ({{ num_contigs }} Contigs)</h3>
            <table class="overview-table">
                <tr><th>Contig</th><th>Inferred Type</th><th>Length</th><th>Circularity</th><th>Evaluation</th></tr>
                {% for c in contig_summaries %}
                <tr>
                    <td><strong>{{ c.id }}</strong></td>
                    <td>{{ c.type }}</td>
                    <td>{{ "%.1f"|format(c.length_kb) }} kb</td>
                    <td><span class="badge-{{ c.j_stat|lower }}">{{ c.j_stat }}</span></td>
                    <td style="color: #555; font-weight: 500;">{{ c.evaluation }}</td>
                </tr>
                {% endfor %}
            </table>
        </div>
        
        <!-- Confidence Profile -->
        <div class="grid">
            <div class="card" style="border-left-color: {% if j_status == 'PASS' %}#2ecc71{% else %}#e74c3c{% endif %};">
                <h3>1. Structural Integrity</h3>
                <div class="value">{{ primary_junc_score }} [{{ j_status }}]</div>
                <div class="desc">Measures physical continuity via read spanning. A passing score confirms the graph is correctly circularized and unbroken.</div>
            </div>
            <div class="card" style="border-left-color: {% if p_status in ['ACCEPTABLE'] %}#2ecc71{% elif p_status in ['WARNING', 'ATYPICAL'] %}#f1c40f{% else %}#e74c3c{% endif %};">
                <h3>2. Replication Structure</h3>
                <div class="value">{{ primary_symmetry }} [{{ p_status }}]</div>
                <div class="desc">Measures bidirectional theta-replication geometry. Evaluates if the assembly represents a biologically sound chromosome.</div>
            </div>
            <div class="card" style="border-left-color: #9b59b6;">
                <h3>3. Base-Level Fidelity</h3>
                <div class="value">{{ motif_count }} Systematic Error(s)</div>
                <div class="desc">Bypasses FASTQ Q-scores by calculating empirical read-to-assembly concordance. Flags specific motifs failed by polishers.</div>
            </div>
        </div>

        <h2>
            Circularity
            <div class="tooltip">?
                <span class="tooltiptext"><strong>Why this matters:</strong> Calculates the ratio of reads that continuously span an artificial assembly junction versus reads that clip/break at that exact coordinate. A score near 0.0 indicates a linear fragment or misassembly, as no physical DNA molecule exists to bridge the gap.</span>
            </div>
            <a href="https://sbam.readthedocs.io/en/latest/Circularity/" target="_blank" class="docs-link">Read Methodology &rarr;</a>
        </h2>
        
        <table>
            <tr><th>Contig</th><th>Length (kb)</th><th>Avg Depth</th><th>Spanning / Broken</th><th>Junction Score</th><th>Status</th></tr>
            {% for contig, metrics in junction.items() %}
            <tr>
                <td>{{ contig }}</td>
                <td>{{ "%.1f"|format(metrics.length / 1000) }}</td>
                <td>{{ metrics.avg_depth }}x</td>
                <td>{{ metrics.spanning_reads }} / {{ metrics.broken_reads }}</td>
                <td>{{ "%.2f"|format(metrics.spanning_score) }}</td>
                <td><span class="badge-{{ metrics.status|lower }}">{{ metrics.status }}</span></td>
            </tr>
            {% endfor %}
        </table>

        <h2>
            Replication structure
            <div class="tooltip">?
                <span class="tooltiptext"><strong>Why this matters:</strong> True bacterial chromosomes replicate bidirectionally, generating a highly symmetrical GC Skew signature pointing to the <i>oriC</i> and <i>ter</i>. SBAM flags contigs that lack this symmetry, highlighting potential misassemblies, chimeras, or unresolvable repeats that standard QC tools miss.</span>
            </div>
            <a href="https://sbam.readthedocs.io/en/latest/Replication-structure/" target="_blank" class="docs-link">Read Methodology &rarr;</a>
        </h2>
        
        <table>
            <tr><th>Contig</th><th>Symmetry</th><th>oriC Pos</th><th>ter Pos</th><th>Plausibility</th></tr>
            {% for contig, metrics in physics.items() %}
            <tr>
                <td>{{ contig }}</td>
                <td>{{ metrics.symmetry }}°</td>
                <td>~{{ metrics.oric_pos }}</td>
                <td>~{{ metrics.ter_pos }}</td>
                <td><span class="badge-{{ metrics.viability|lower }}">{{ metrics.viability }}</span></td>
            </tr>
            {% endfor %}
            {% if not physics %}
            <tr><td colspan="5" style="text-align: center;">No contigs >300kb available for replication architecture analysis.</td></tr>
            {% endif %}
        </table>

        <!-- Interactive Plotly Render Area -->
        {% for contig, metrics in physics.items() %}
        <div class="plot-container">
            <h3 style="text-align: center; color: #2c3e50; margin-bottom: 5px;">Genome Architecture Map: Contig {{ contig }}</h3>
            <div id="plot_{{ contig }}" style="width:100%; max-width:800px; height:600px; margin:0 auto;"></div>

                       <div style="text-align: center; margin-top: 10px; margin-bottom: 20px; font-size: 13px; color: #7f8c8d; display: flex; justify-content: center; gap: 20px; flex-wrap: wrap;">
                <div style="display: flex; align-items: center;">
                    <span style="display: inline-block; width: 20px; height: 6px; background-color: #bdc3c7; margin-right: 8px; border-radius: 3px;"></span>
                    Genome Track
                    <div class="tooltip">?
                        <span class="tooltiptext">The circle represents the total length of the contig.</span>
                    </div>
                </div>
                <div style="display: flex; align-items: center;">
                    <span style="display: inline-block; width: 20px; height: 6px; background-color: #3498db; margin-right: 8px; border-radius: 3px; opacity: 0.7;"></span>
                    GC Skew
                    <div class="tooltip">?
                        <span class="tooltiptext">Tracks the cumulative Guanine vs. Cytosine bias. In a true chromosome, bidirectional replication causes the leading strand to enrich with Guanine, creating a distinctive smooth wave.</span>
                    </div>
                </div>
                <div style="display: flex; align-items: center;">
                    <span style="display: inline-block; width: 10px; height: 10px; background-color: #2ecc71; border-radius: 50%; margin-right: 4px;"></span>
                    <span style="display: inline-block; width: 10px; height: 10px; background-color: #e74c3c; border-radius: 50%; margin-right: 8px;"></span>
                    Key Loci (oriC & ter)
                    <div class="tooltip">?
                        <span class="tooltiptext">The predicted Origin of Replication and Terminus. On a perfectly assembled chromosome, these should be physically opposite each other (~180° apart).</span>
                    </div>
                </div>
                <div style="display: flex; align-items: center;">
                    <span style="display: inline-block; width: 20px; height: 2px; border-top: 2px dotted #e74c3c; margin-right: 8px;"></span>
                    Skew = 0
                    <div class="tooltip">?
                        <span class="tooltiptext">The mathematical baseline (zero-crossing) of the cumulative GC skew calculation.</span>
                    </div>
                </div>
            </div>
            
            <script>
                var raw_r = {{ metrics.plot_r | tojson }};
                var theta = {{ metrics.plot_theta | tojson }};
                var oric_pos = {{ metrics.oric_pos }};
                var ter_pos = {{ metrics.ter_pos }};
                var seq_len = {{ metrics.length }};

                var min_r = Math.min(...raw_r);
                var max_r = Math.max(...raw_r);
                var range = (max_r - min_r) || 1;
                var baseline = 2.0;
                
                var norm_r = raw_r.map(r => baseline + ((r - min_r) / range - 0.5) * 1.6);
                var custom_data = raw_r.map((val, i) => {
                    var bp = Math.round((theta[i] / 360) * seq_len);
                    return [bp.toLocaleString(), Math.round(val).toLocaleString()];
                });

                var contig_trace = { type: 'scatterpolar', r: Array(theta.length).fill(baseline), theta: theta, mode: 'lines', line: {color: '#bdc3c7', width: 6}, hoverinfo: 'none', name: 'Genome Track' };

                var skew_trace = {
                    type: 'scatterpolar', r: norm_r, theta: theta, mode: 'lines', line: {color: '#3498db', width: 2}, fill: 'tonext', fillcolor: 'rgba(52, 152, 219, 0.3)',
                    customdata: custom_data, hovertemplate: '<b>Pos:</b> %{customdata[0]} bp<br><b>Cum. Skew:</b> %{customdata[1]}<extra></extra>',
                    name: 'GC Skew (Min: ' + Math.round(min_r).toLocaleString() + ', Max: ' + Math.round(max_r).toLocaleString() + ')'
                };

                var oric_angle = (oric_pos / seq_len) * 360;
                var ter_angle = (ter_pos / seq_len) * 360;
                var marker_trace = {
                    type: 'scatterpolar', r: [baseline, baseline], theta: [oric_angle, ter_angle], mode: 'markers+text',
                    marker: { color: ['#2ecc71', '#e74c3c'], size: 14, line: {color: 'white', width: 2} },
                    text: ['<b>oriC</b>', '<b>ter</b>'], textposition: ['top right', 'bottom left'], textfont: {size: 16, color: '#2c3e50'},
                    hoverinfo: 'text', hovertext: [`oriC: ~${oric_pos.toLocaleString()} bp`, `ter: ~${ter_pos.toLocaleString()} bp`], name: 'Key Loci'
                };

                var r_tickvals = [baseline - 0.8, baseline, baseline + 0.8];
                var mid_val = (min_r + max_r) / 2;
                var r_ticktext = [ Math.round(min_r).toLocaleString(), Math.round(mid_val).toLocaleString(), Math.round(max_r).toLocaleString() ];
                var plot_data = [contig_trace, skew_trace, marker_trace];

                if (min_r < 0 && max_r > 0) {
                    var zero_norm = baseline + ((0 - min_r) / range - 0.5) * 1.6;
                    r_tickvals.push(zero_norm); r_ticktext.push("0");
                    plot_data.push({ type: 'scatterpolar', r: Array(theta.length).fill(zero_norm), theta: theta, mode: 'lines', line: {color: '#e74c3c', width: 1, dash: 'dot'}, hoverinfo: 'none', name: 'Skew = 0' });
                }

                var tick_vals = []; var tick_text = [];
                for (var i = 0; i < 8; i++) {
                    tick_vals.push((i * 360) / 8);
                    tick_text.push(Math.round((i * seq_len) / 8 / 1000).toLocaleString() + ' kb');
                }

                var layout = {
                    polar: {
                        radialaxis: { visible: true, range: [0, 4], showticklabels: true, tickvals: r_tickvals, ticktext: r_ticktext, angle: 90, tickangle: 0, tickfont: {size: 11, color: '#e67e22', weight: 'bold'}, showgrid: true, gridcolor: 'rgba(0,0,0,0.05)', gridwidth: 1 },
                        angularaxis: { direction: "clockwise", visible: true, tickmode: 'array', tickvals: tick_vals, ticktext: tick_text, showgrid: true, gridcolor: 'rgba(0,0,0,0.1)', linecolor: 'rgba(0,0,0,0)' }
                    },
                    showlegend: false, margin: { t: 60, b: 20, l: 60, r: 60 }
                };

                Plotly.newPlot('plot_{{ contig }}', plot_data, layout, {responsive: true});
            </script>
        </div>
        {% endfor %}

        <h2 style="margin-top: 50px;">
            Base-Level Analysis
            <div class="tooltip">?
                <span class="tooltiptext"><strong>The Limitation of Q-Scores:</strong> Consensus Q-scores generated by polishers are algorithmic heuristics, not physical measurements. SBAM maps the original physical sequencing reads back to the FASTA to empirically measure true sequence concordance and flag systematic motif dropouts caused by basecaller errors.</span>
            </div>
            <a href="https://sbam.readthedocs.io/en/latest/Base-Level-Analysis/" target="_blank" class="docs-link">Read Methodology &rarr;</a>
        </h2>
        
        <!-- Dynamic Explanation Box -->
        {% if masked_pct > 5.0 %}
        <div class="explanation-box-warning">
            <h4>Excluded Structural Regions ({{ masked_bases }} bp / {{ "%.2f"|format(masked_pct) }}%)</h4>
            <p style="margin-top: 10px;">A fraction of masked regions &gt;5% indicates a fragmented assembly or unresolved repeats, lowering confidence in the genome's completeness. <strong>{{ "%.2f"|format(evaluated_pct) }}%</strong> of the primary chromosome was successfully evaluated for motif fidelity.</p>
        </div>
        {% else %}
        <div class="explanation-box-success">
            <h4>Excluded Structural Regions ({{ masked_bases }} bp / {{ "%.2f"|format(masked_pct) }}%)</h4>
            <p style="margin-top: 10px;">A fraction of masked regions &lt;5% indicates a highly intact assembly, increasing confidence in the genome's completeness. <strong>{{ "%.2f"|format(evaluated_pct) }}%</strong> of the primary chromosome was successfully evaluated for motif fidelity.</p>
        </div>
        {% endif %}
        
        {% if masked_regions %}
        <h3 style="font-size: 14px; color: #34495e; margin-top: 20px;">
        Largest Excluded Regions ( 50 bp)
        <div class="tooltip">?
            <span class="tooltiptext">
            <strong>High Depth:</strong> Indicates unresolved collapsed repeats<br>
            <strong>Low Depth/Discordance:</strong> Indicates gaps, misassemblies
            </span>
        </div>        
        </h3>
        <table style="margin-bottom: 30px; font-size: 13px;">
            <tr><th>Region Start</th><th>Region End</th><th>Length (bp)</th><th>Primary Reason</th></tr>
            {% for region in masked_regions %}
            <tr>
                <td>{{ region.start }}</td>
                <td>{{ region.end }}</td>
                <td><strong>{{ region.length }}</strong></td>
                <td style="color: #e67e22;">{{ region.reason }}</td>
            </tr>
            {% endfor %}
        </table>
        {% endif %}
        
        <h3 style="font-size: 14px; color: #34495e; margin-top: 20px;">Systematic Motif Errors (Evaluated Regions)</h3>
        <table>
            <tr><th>k-mer Size</th><th>Motif</th><th>Discordant Count / Total</th><th>Error Rate</th><th>Fold Increase</th></tr>
            {% for motif in motifs %}
            <tr>
                <td>{{ motif.kmer_size }}</td>
                <td><strong>{{ motif.motif }}</strong></td>
                <td>{{ motif.discordant_count }} / {{ motif.occurrences }}</td>
                <td>{{ "%.1f"|format(motif.error_rate * 100) }}%</td>
                <td>{{ "%.1f"|format(motif.fold_increase) }}x</td>
            </tr>
            {% endfor %}
            {% if not motifs %}
            <tr><td colspan="5" style="text-align: center; color: #27ae60; font-weight: bold; padding: 20px;">No systematic motif errors detected in evaluated regions. Basecalling is highly concordant!</td></tr>
            {% endif %}
        </table>
    </div>
</body>
</html>
"""

class DashboardBuilder:
    def __init__(self, outdir, assembly_path):
        self.outdir = outdir
        self.assembly_name = os.path.basename(assembly_path)
        
    def generate_report(self, junction_metrics, physics_metrics, motif_results, masked_bases, masked_pct, masked_regions=None):
        print(" > [Report] Compiling HTML dashboard...")
        
        # identify chromosome
        primary_contig_id = None
        if junction_metrics:
            primary_contig_id = max(junction_metrics.keys(), key=lambda k: junction_metrics[k]['length'])
        
        # calculate confidence
        primary_junc_score = "N/A"
        j_status = "UNKNOWN"
        if primary_contig_id and primary_contig_id in junction_metrics:
            j_status = junction_metrics[primary_contig_id]['status']
            primary_junc_score = f"{junction_metrics[primary_contig_id]['spanning_score']:.2f}"
            
        primary_symmetry = "N/A"
        p_status = "N/A"
        if primary_contig_id and primary_contig_id in physics_metrics:
            p_status = physics_metrics[primary_contig_id]['viability']
            primary_symmetry = f"{physics_metrics[primary_contig_id]['symmetry']}°"
            
        motif_count = len(motif_results) if motif_results else 0
        evaluated_pct = 100.0 - masked_pct
        
        # Contig summary
        num_contigs = len(junction_metrics) if junction_metrics else 0
        contig_summaries = []
        if junction_metrics:
            for cid, j_metrics in junction_metrics.items():
                length_kb = j_metrics['length'] / 1000
                c_j_stat = j_metrics['status']
                c_type = "Chromosome" if j_metrics['length'] >= 1000000 else "Plasmid/Fragment"
                c_p_stat = physics_metrics.get(cid, {}).get('viability', 'N/A')
                
                # Dynamic biological evaluation
                if c_j_stat == "PASS":
                    if c_p_stat in ["ACCEPTABLE"]:
                        evaluation = "Circular"
                    elif c_p_stat in ["ATYPICAL", "WARNING"]:
                        evaluation = "Structurally Intact, Atypical replication structure"
                    else:
                        evaluation = "Structurally Intact Plasmid"
                else:
                    evaluation = "Linear Fragment or Misassembly"
                    
                contig_summaries.append({
                    "id": cid,
                    "type": c_type,
                    "length_kb": length_kb,
                    "j_stat": c_j_stat,
                    "evaluation": evaluation
                })
        
        # output result
        #result = "FAIL"
        result_color = "#c0392b"
        result_msg = "Assembly fails primary structural spanning thresholds. High probability of misassembly, " \
        "linear fragmentation, or unresolved repeats at contig boundaries."
        
        if j_status == "PASS":
            if p_status in ["ACCEPTABLE"]:
                #result = "ACCEPT"
                result_color = "#27ae60"
                result_msg = "Strong evidence of structural circularity and expected biological replication architecture. " \
                "High confidence in chromosomal integrity."
            elif p_status in ["ATYPICAL", "WARNING"]:
                #result = "REVIEW"
                result_color = "#f39c12"
                result_msg = "Structurally intact, but the replication architecture deviates significantly from standard bidirectional " \
                "theta-replication. May indicate a chimera, major inversion, or highly atypical biology."
            else:
                #result = "ACCEPT WITH CAUTION"
                result_color = "#f39c12"
                result_msg = "Structurally intact, but chromosomal size is too small to perform replication analysis. " \
                "Likely a plasmid assembly."

        template = Template(HTML_TEMPLATE)
        html_content = template.render(
            assembly_name=self.assembly_name,
            date=datetime.now().strftime("%Y-%m-%d %H:%M"),
            #result=result,
            result_color=result_color,
            result_msg=result_msg,
            num_contigs=num_contigs,
            contig_summaries=contig_summaries,
            j_status=j_status,
            p_status=p_status,
            primary_junc_score=primary_junc_score,
            primary_symmetry=primary_symmetry,
            motif_count=motif_count,
            junction=junction_metrics,
            physics=physics_metrics,
            motifs=motif_results,
            masked_bases=masked_bases,
            masked_pct=masked_pct,
            evaluated_pct=evaluated_pct,
            masked_regions=masked_regions
        )
        
        report_path = os.path.join(self.outdir, f"{self.assembly_name}_sbam_report.html")
        with open(report_path, "w") as f:
            f.write(html_content)
            
        print(f" > [Report] Dashboard saved to: {report_path}")
        return report_path