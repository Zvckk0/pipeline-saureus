import pandas as pd
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

print("\n" + "="*70)
print("ARBITRE AMR : CONFRONTATION ABRICATE vs AMRFINDER PLUS")
print("="*70 + "\n")

try:
    # 1. Extraction et filtrage ABRicate (Base CARD)
    df_ab = pd.read_csv("resultats_card.tsv", sep='\t')
    df_ab_filtre = df_ab[(df_ab['%COVERAGE'] >= 90.0) & (df_ab['%IDENTITY'] >= 90.0)]
    genes_abricate = set(df_ab_filtre['GENE'].dropna().tolist())
    
    # 2. Extraction et filtrage AMRFinder Plus
    df_amr = pd.read_csv("resultats_amrfinder.tsv", sep='\t')
    # Les noms de colonnes corrigés pour AMRFinder
    df_amr_filtre = df_amr[(df_amr['% Coverage of reference'] >= 90.0) & (df_amr['% Identity to reference'] >= 90.0)]
    genes_amrfinder = set(df_amr_filtre['Element symbol'].dropna().tolist())

    # 3. Synthèse des résultats
    commun = genes_abricate.intersection(genes_amrfinder)
    
    print(f"Gènes identifiés par ABRicate (CARD) : {len(genes_abricate)}")
    print(f"Gènes identifiés par AMRFinder Plus : {len(genes_amrfinder)}")
    print(f"\nGènes validés par les DEUX algorithmes (Haute Confiance) : {len(commun)}")
    
    if commun:
        print(" -> Liste des gènes confirmés :", ", ".join(commun))

    # 4. Création du Diagramme de Venn
    plt.figure(figsize=(8, 6))
    venn2([genes_abricate, genes_amrfinder], ('ABRicate (CARD)', 'AMRFinder Plus'))
    plt.title("Croisement des gènes de résistance (S. aureus)")
    
    # Sauvegarde de l'image
    fichier_sortie = "diagramme_venn_amr.png"
    plt.savefig(fichier_sortie, dpi=300, bbox_inches='tight')
    print(f"\n[SUCCÈS] Le diagramme a été généré et sauvegardé sous : {fichier_sortie}\n")

except FileNotFoundError as e:
    print(f"Erreur de lecture : {e}. Vérifie que les fichiers .tsv sont bien dans ce dossier.")
except Exception as e:
    print(f"Une erreur inattendue s'est produite : {e}")
