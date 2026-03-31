# ------------------------------------------------------------
# Figure 1 — Patient selection flowchart
# Same structure, with detailed counts added inside combined boxes
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
E4 <- N3 - N8   # Combined eligibility filters
E5 <- N8 - N10  # Missing/unknown covariates or staging

# Detailed counts within combined exclusion steps
E4a <- 10     # Age <18 years
E4b <- 20081  # Not first malignant primary
E4c <- 3443   # >1 malignant tumor
E4d <- 13514  # Missing/unknown/0 survival months
E4e <- 10281  # Unknown metastasis indicators

E5a <- 44175  # Missing/unknown covariates or C34.9
E5b <- 2124   # TX or NX stage

g <- grViz(sprintf('
digraph flow {
  graph [rankdir=TB, nodesep=0.40, ranksep=0.65, fontsize=12]
  node  [shape=rectangle, fontname=Helvetica, fontsize=12,
         style="rounded,filled", fillcolor="white", color="black",
         penwidth=1.2, width=3]
  edge  [fontname=Helvetica, fontsize=11, color="black", arrowsize=0.9]

  // Main vertical flow
  A [label="Lung cancer (SEER)\\nN = %s"];
  B [label="Diagnosis years\\n2004–2016\\nN = %s"];
  C [label="Non-small cell lung cancer\\nN = %s"];
  D [label="Metastatic at diagnosis (M1)\\nN = %s"];
  E [label="Eligibility filters applied\\nN = %s"];
  F [label="Final metastatic NSCLC cohort\\nN = %s", penwidth=2.0];

  // Right-hand exclusion boxes
  node [shape=rectangle, style="rounded,filled",
        fillcolor="white", color="black", fontsize=10.5, width=4.2]

  Bx [label="Excluded:\\nOutside 2004–2016\\n(n = %s)"];
  Cx [label="Excluded:\\nSmall-cell lung cancer\\n(n = %s)"];
  Dx [label="Excluded:\\nM0, MX, or unknown M stage\\n(n = %s)"];
  Ex [label="Excluded:\\nAge <18 years (n = %s); prior malignant tumor (n = %s);\\n>1 malignant tumor (n = %s); missing/0 survival months (n = %s);\\nunknown metastasis indicators (n = %s)\\nTotal excluded at this step: n = %s"];
  Fx [label="Excluded:\\nMissing/unknown covariates or C34.9-Lung, NOS (n = %s);\\nTX or NX stage (n = %s)\\nTotal excluded at this step: n = %s"];

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
fmt(E1), fmt(E2), fmt(E3),
fmt(E4a), fmt(E4b), fmt(E4c), fmt(E4d), fmt(E4e), fmt(E4),
fmt(E5a), fmt(E5b), fmt(E5)
))

# Preview in R
g

# Export SVG + PDF
svg_txt <- export_svg(g)
cat(svg_txt, file = "Figure_flowchart_wrap_revised.svg")

rsvg_pdf("Figure_flowchart_wrap_revised.svg",
         file  = "Figure_flowchart_wrap_revised_doublecol.pdf",
         width = 7.2,
         height = 11)

cat("Revised flowchart saved successfully.\n")