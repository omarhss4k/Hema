# Diagnostic — structure réelle de PK_singe_FGFR2.xlsx
library(readxl)

raw <- suppressMessages(
  read_xlsx("PK_singe_FGFR2.xlsx",
            sheet = "Sheet1",
            col_names = FALSE,
            .name_repair = "minimal")
)

cat("=== Dimensions :", nrow(raw), "lignes ×", ncol(raw), "colonnes ===\n\n")

cat("=== 30 premières lignes (colonnes 1-6) ===\n")
print(as.data.frame(raw[1:min(30, nrow(raw)), 1:min(6, ncol(raw))]))

cat("\n=== Contenu unique de la colonne 1 ===\n")
print(unique(as.character(raw[[1]])))

cat("\n=== Lignes contenant 'Time' ou 'time' (col 1) ===\n")
idx <- grep("(?i)time", as.character(raw[[1]]))
cat("Indices :", idx, "\n")
if (length(idx) > 0) {
  for (i in idx) {
    cat("Ligne", i, ":", paste(as.character(unlist(raw[i, 1:min(8, ncol(raw))])), collapse = " | "), "\n")
  }
}

cat("\n=== Noms des feuilles disponibles ===\n")
print(excel_sheets("PK_singe_FGFR2.xlsx"))
