# ------------------------------------------------------------
# Figure 1 — Patient selection flowchart (B&W, journal-ready)
# ------------------------------------------------------------
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

# Counts
N0  <- 1055564  # Lung cancer (SEER)
N1  <- 628425   # Diagnosis years 2004–2016
N2  <- 553763   # NSCLC
N3  <- 112882   # M1 at diagnosis
N8  <- 65553    # Eligibility filters applied
N10 <- 19254    # Final cohort

fmt <- function(x) formatC(x, format = "d", big.mark = ",")

E1 <- N0 - N1   # Outside 2004–2016
E2 <- N1 - N2   # SCLC
E3 <- N2 - N3   # M0 / MX / unknown M
E4 <- N3 - N8   # Age <18, not first primary, >1 tumor, 0 months, unknown mets
E5 <- N8 - N10  # Missing/unknown covariates or staging

g <- grViz(sprintf('
digraph flow {
  graph [rankdir=TB, nodesep=0.40, ranksep=0.65, fontsize=12]
  node  [shape=rectangle, fontname=Helvetica, fontsize=12,
         style="rounded,filled", fillcolor="white", color="black",
         penwidth=1.2, width=3]   // narrower width so text fits better
  edge  [fontname=Helvetica, fontsize=11, color="black", arrowsize=0.9]

  // Main vertical flow
  A [label="Lung cancer (SEER)\\nN = %s"];
  B [label="Diagnosis years\\n2004–2016\\nN = %s"];
  C [label="Non-small cell lung cancer\\nN = %s"];
  D [label="Metastatic at diagnosis (M1)\\nN = %s"];
  E [label="Eligibility filters applied\\nN = %s"];
  F [label="Final metastatic NSCLC cohort\\nN = %s", penwidth=2.0]

  // Right-hand exclusion boxes
  node [shape=rectangle, style="rounded,filled",
        fillcolor="white", color="black", fontsize=11, width=3]

  Bx [label="Excluded:\\nOutside 2004–2016\\n(n = %s)"];
  Cx [label="Excluded:\\nSmall-cell lung cancer\\n(n = %s)"];
  Dx [label="Excluded:\\nM0, MX, or unknown M stage\\n(n = %s)"];
  Ex [label="Excluded:\\nAge < 18 years; prior malignant tumor;\\n>1 malignant tumor; 0 survival months; or\\nunknown metastasis indicators\\n(n = %s)"];
  Fx [label="Excluded:\\nMissing or unknown covariates\\nor staging information\\n(n = %s)"];

  // Flow arrows
  A -> B -> C -> D -> E -> F;

  // Same-rank alignment
  {rank=same; B; Bx}
  {rank=same; C; Cx}
  {rank=same; D; Dx}
  {rank=same; E; Ex}
  {rank=same; F; Fx}

  // Dashed connectors
  edge [style=dashed, arrowsize=0.6, color="black", constraint=false]
  B -> Bx [dir=none]
  C -> Cx [dir=none]
  D -> Dx [dir=none]
  E -> Ex [dir=none]
  F -> Fx [dir=none]

  // Invisible spine for right column
  edge [style=invis, constraint=true]
  Bx -> Cx -> Dx -> Ex -> Fx;
}
',
fmt(N0), fmt(N1), fmt(N2), fmt(N3), fmt(N8), fmt(N10),
fmt(E1), fmt(E2), fmt(E3), fmt(E4), fmt(E5)
))

# Preview in R
g

# Export SVG + PDF
svg_txt <- export_svg(g)
cat(svg_txt, file = "Figure1_flowchart_wrap.svg")

rsvg_pdf("Figure1_flowchart_wrap.svg",
         file  = "Figure1_flowchart_wrap_doublecol.pdf",
         width = 7.2,   # double-column width
         height = 11)   # enough vertical space
