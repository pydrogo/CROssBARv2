import pandas as pd
import numpy as np


def node_creator():
    """
    This function takes random identifier and name(if available) samples
    from datasets and creates dataframes from these samples. After,
    these dataframes were concatenated into single dataset.

    Sizes: gene_list_df -  300 protein_list_df - 300 disease_list_df -
    200 pathway_list_df - 200 drug_list_df - 100 compound_list_df - 100
    phenotype_list_df - 100 function_list_df - 100 tissue_list_df - 200
    patient_list_df - 200 domain_list_df - 100 side_effect_list_df - 100
    location_list_df - 100 organism_list_df - 1

    TOTAL: 1959 rows × 4 columns
    """
    # Gene and Protein ids
    gene_name_and_protein_acc_df = pd.read_csv(
        r"data/datasets/gene_name_and_protein_acc.csv"
    )
    sampled_gene_and_protein = gene_name_and_protein_acc_df.sample(
        n=300, random_state=49
    ).reset_index(drop=True)

    gene_list = []
    protein_list = []
    for idx, row in sampled_gene_and_protein.iterrows():
        gene_list.append(("Gene", "UniprotKB", row["gene"], np.nan))
        protein_list.append(("Protein", "UniprotKB", row["acc"], np.nan))

    gene_list_df = pd.DataFrame(
        gene_list, columns=["Type", "Source Database", "ID", "Name"]
    )
    protein_list_df = pd.DataFrame(
        protein_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # Efo disease ids and names
    efo_diseases_df = pd.read_csv(r"data/datasets/efo_disease_terms.csv")
    sampled_efo = efo_diseases_df.sample(n=100, random_state=49).reset_index(
        drop=True
    )

    efo_diseases_list = []
    for idx, row in sampled_efo.iterrows():
        efo_diseases_list.append(
            ("Disease", "EFO", row["obo_id"], row["label"])
        )

    efo_diseases_list_df = pd.DataFrame(
        efo_diseases_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # KEGG pathway and disease ids and names
    kegg_pathway_and_disease_df = pd.read_csv(
        r"data/datasets/kegg_pathway_and_disease.csv"
    )
    kegg_pathway_and_disease_df = kegg_pathway_and_disease_df.drop_duplicates(
        subset=["kegg_pathwayname", "kegg_diseasename"]
    )
    sampled_kegg_pathway_and_disease = kegg_pathway_and_disease_df.sample(
        n=100, random_state=49
    ).reset_index(drop=True)

    kegg_pathway_list = []
    kegg_disease_list = []
    for idx, row in sampled_kegg_pathway_and_disease.iterrows():
        kegg_pathway_list.append(
            ("Pathway", "KEGG", row["kegg_pathwayid"], row["kegg_pathwayname"])
        )
        kegg_disease_list.append(
            ("Disease", "KEGG", row["kegg_diseaseid"], row["kegg_diseasename"])
        )

    kegg_pathway_list_df = pd.DataFrame(
        kegg_pathway_list, columns=["Type", "Source Database", "ID", "Name"]
    )
    kegg_disease_list_df = pd.DataFrame(
        kegg_disease_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # Reactome pathway ids and names
    reactome_pathway_df = pd.read_csv(r"data/datasets/reactome_pathway.csv")
    sampled_reactome_pathway = reactome_pathway_df.sample(
        n=100, random_state=49
    ).reset_index(drop=True)

    reactome_pathway_list = []
    for idx, row in sampled_reactome_pathway.iterrows():
        reactome_pathway_list.append(
            ("Pathway", "Reactome", row["id"], row["pathwayName"])
        )

    reactome_pathway_list_df = pd.DataFrame(
        reactome_pathway_list,
        columns=["Type", "Source Database", "ID", "Name"],
    )

    # concat disease (efo and kegg) and pathway (reactome and kegg) dataframes
    disease_list_df = pd.concat(
        [efo_diseases_list_df, kegg_disease_list_df]
    ).reset_index(drop=True)
    pathway_list_df = pd.concat(
        [reactome_pathway_list_df, kegg_pathway_list_df]
    ).reset_index(drop=True)

    # drug ids and names
    drug_df = pd.read_csv(r"data/datasets/drug_terms.csv")
    sampled_drug_df = drug_df.sample(n=100, random_state=49).reset_index(
        drop=True
    )

    drug_list = []
    for idx, row in sampled_drug_df.iterrows():
        drug_list.append(("Drug", "DrugBank", row["identifier"], row["name"]))

    drug_list_df = pd.DataFrame(
        drug_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # compound ids
    compound_df = pd.read_csv(r"data/datasets/compound_terms.csv")
    sampled_compound_df = compound_df.sample(
        n=100, random_state=49
    ).reset_index(drop=True)

    compound_list = []
    for idx, row in sampled_compound_df.iterrows():
        compound_list.append(
            ("Compound", "ChEMBL", row["molecule_chembl_id"], np.nan)
        )

    compound_list_df = pd.DataFrame(
        compound_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # Phenotype ids and names
    phenotype_df = pd.read_csv(r"data/datasets/phenotype_terms.csv")
    sampled_phenotype_df = phenotype_df.sample(
        n=100, random_state=49
    ).reset_index(drop=True)

    phenotype_list = []
    for idx, row in sampled_phenotype_df.iterrows():
        phenotype_list.append(
            ("Phenotype", "HPO", row["hpo_id"], row["term_name"])
        )

    phenotype_list_df = pd.DataFrame(
        phenotype_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # go function ids and names
    go_function_df = pd.read_csv(r"data/datasets/go_function_terms.csv")
    go_function_df = go_function_df.drop_duplicates(
        subset=["UniProt", "GO Term go_id"]
    )
    go_function_df.dropna(axis=0, inplace=True)
    sampled_go_function_df = go_function_df.sample(
        n=100, random_state=49
    ).reset_index(drop=True)

    function_list = []
    for idx, row in sampled_go_function_df.iterrows():
        function_list.append(
            ("Function", "GO", row["GO Term go_id"], row["GO Term"])
        )

    function_list_df = pd.DataFrame(
        function_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # Cell/Tissue ids and names
    tissue_terms_df = pd.read_csv(
        r"data/datasets/tissue_terms.csv"
    ).reset_index()
    tissue_terms_df.drop(labels=["level_0"], axis=1, inplace=True)
    tissue_terms_df.columns = [
        "Name",
        "Cell Model Passports",
        "COSMIC ID",
        "TCGA Classification",
        "Tissue",
        "Tissue sub-type",
        "COUNT1",
        "COUNT2",
        "COUNT3",
    ]
    tissue_terms_df = tissue_terms_df.drop_duplicates(
        subset="Cell Model Passports"
    )
    sampled_tissue_terms_df = tissue_terms_df.sample(
        n=200, random_state=49
    ).reset_index(drop=True)

    tissue_list = []
    for idx, row in sampled_tissue_terms_df.iterrows():
        tissue_list.append(
            ("Cell/Tissue", "GDSC", row["Cell Model Passports"], row["Tissue"])
        )

    tissue_list_df = pd.DataFrame(
        tissue_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # patient ids and names
    frames = []
    for n in range(5):
        patient_df = pd.read_csv(
            rf"data/datasets/patient_terms/repository-cases-table.2022-07-09 ({n}).tsv",
            sep="\t",
        )
        frames.append(patient_df)
    patient_df = pd.concat(frames)
    sampled_patient_df = patient_df.sample(n=200, random_state=49).reset_index(
        drop=True
    )

    patient_list = []
    for idx, row in sampled_patient_df.iterrows():
        patient_list.append(
            ("Patient", "TCGA", row["Case ID"], row["Primary Site"])
        )

    patient_list_df = pd.DataFrame(
        patient_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # domain ids and names
    domain_df = pd.read_csv(r"data/datasets/domain_terms.tsv", sep="\t")
    domain_df.drop(
        labels=["Integrated Signatures", "GO Terms"], axis=1, inplace=True
    )
    domain_df.dropna(axis=0, inplace=True)
    domain_df = domain_df[domain_df["Type"] == "domain"]
    sampled_domain_df = domain_df.sample(n=100, random_state=49).reset_index(
        drop=True
    )

    domain_list = []
    for idx, row in sampled_domain_df.iterrows():
        domain_list.append(
            ("Domain", "InterPro", row["Integrated Into"], row["Name"])
        )

    domain_list_df = pd.DataFrame(
        domain_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # side effect ids and names
    side_effect_df = pd.read_excel(r"data/datasets/Drug_ADR.xlsx")
    sampled_side_effect_df = side_effect_df.sample(
        n=100, random_state=49
    ).reset_index(drop=True)

    side_effect_list = []
    for idx, row in sampled_side_effect_df.iterrows():
        side_effect_list.append(
            ("Side Effect", "ADReCS", row["ADR_ID"], row["ADR_TERM"])
        )

    side_effect_list_df = pd.DataFrame(
        side_effect_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # location ids and names
    location_df = pd.read_csv(r"data/datasets/location_terms.tsv", sep="\t")
    sampled_location_df = location_df.sample(
        n=100, random_state=49
    ).reset_index(drop=True)

    location_list = []
    for idx, row in sampled_location_df.iterrows():
        location_list.append(
            ("Location", "GOA", row["ECO ID"], row["GO NAME"])
        )

    location_list_df = pd.DataFrame(
        location_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # organism id and name
    organism_list = [("Organism", "UniprotKB", "9606", "Homo sapiens")]
    organism_list_df = pd.DataFrame(
        organism_list, columns=["Type", "Source Database", "ID", "Name"]
    )

    # Finally concatenate them
    all_nodes_list = [
        gene_list_df,
        protein_list_df,
        disease_list_df,
        pathway_list_df,
        drug_list_df,
        compound_list_df,
        phenotype_list_df,
        function_list_df,
        tissue_list_df,
        patient_list_df,
        domain_list_df,
        side_effect_list_df,
        location_list_df,
        organism_list_df,
    ]

    all_nodes_list_df = pd.concat(all_nodes_list).reset_index(drop=True)
    all_nodes_list_df.drop_duplicates(inplace=True)

    # write to csv
    all_nodes_list_df.to_csv(r"data/nodes/nodes.csv", index=False)

    # return final dataframe
    return all_nodes_list_df


def create_relationship(
    df, source, target, label, interaction_number, randomly_choose=True
):
    """
    Creates interaction pairs from nodes dataframe
    args:
        df -> nodes_df
        source -> source term of directed edge (node type in the dataframe)
        target -> target term of directed edge (node type in the dataframe)
        label -> label of edge
        interaction_number[int] -> total interaction number between pairs
        randomly choose[boolean] -> whether randomly select source and target or not
    """
    source_df = df[df["Type"] == source].reset_index(drop=True)
    target_df = df[df["Type"] == target].reset_index(drop=True)

    for count in range(interaction_number):
        if randomly_choose:
            while True:
                source_choice = np.random.choice(
                    source_df.index.values, size=1
                )[0]
                target_choice = np.random.choice(
                    target_df.index.values, size=1
                )[0]
                if source_choice != target_choice:
                    break
        else:
            if target == "Organism":
                source_choice = count
                target_choice = 0
            else:
                source_choice = count
                target_choice = count

        source_id = source_df.iloc[source_choice, :]["ID"]
        target_id = target_df.iloc[target_choice, :]["ID"]

        yield (source, source_id, target, target_id, label)


if __name__ == "__main__":

    nodes_df = node_creator()

    # mapping of interactions
    Interactions_map = [
        ("Protein", "Protein", "Interacts_With", 150, True),
        ("Protein", "Domain", "Has", 85, True),
        ("Protein", "Function", "Enables", 85, True),
        ("Protein", "Location", "Localizes_To", 85, True),
        ("Protein", "Pathway", "Take_Part_In", 175, True),
        ("Gene", "Protein", "Encodes", 300, False),
        ("Gene", "Gene", "Is_Ortholog_To", 100, True),
        ("Gene", "Gene", "Regulates", 100, True),
        ("Gene", "Organism", "Belongs_To", 300, False),
        ("Gene", "Cell/Tissue", "Is_Mutated_In", 85, True),
        ("Gene", "Cell/Tissue", "Is_DEG_In", 85, True),
        ("Gene", "Phenotype", "Is_Associated_With", 85, True),
        ("Gene", "Patient", "Is_Mutated_In", 85, True),
        ("Gene", "Patient", "Is_DEG_In", 85, True),
        ("Gene", "Pathway", "Is_Member_Of", 175, True),
        ("Gene", "Disease", "Is_related_to", 175, True),
        ("Disease", "Protein", "Is_related_to", 175, True),
        ("Drug", "Drug", "Interacts_With", 75, True),
        ("Drug", "Protein", "Targets", 85, True),
        ("Drug", "Side Effect", "Has", 75, True),
        ("Drug", "Cell/Tissue", "Targets", 85, True),
        ("Drug", "Pathway", "Has_Target_In", 85, True),
        ("Compound", "Protein", "Targets", 85, True),
        ("Cell/Tissue", "Disease", "Has", 85, True),
        ("Patient", "Disease", "Has", 100, True),
        ("Disease", "Disease", "Comorbid_With", 100, True),
        ("Disease", "Phenotype", "Is_Associated_With", 75, True),
        ("Disease", "Pathway", "Modulates", 85, True),
        ("Disease", "Drug", "Is_Treated_By", 50, True),
        ("Domain", "Function", "Has", 50, True),
        ("Domain", "Location", "Has", 50, True),
    ]

    all_interactions = []
    for relation in Interactions_map:
        print(relation)
        interactions = list(
            create_relationship(
                df=nodes_df,
                source=relation[0],
                target=relation[1],
                label=relation[2],
                interaction_number=relation[3],
                randomly_choose=relation[4],
            )
        )
    all_interactions.extend(interactions)

    # create dataframe from all_interactions
    edges_df = pd.DataFrame(
        all_interactions,
        columns=[
            "Source Type",
            "Source ID",
            "Target Type",
            "Target ID",
            "Label",
        ],
    )
    edges_df.drop_duplicates(inplace=True)

    # write to csv
    edges_df.to_csv(
        r"data/edges/edges.csv", index=False
    )  # size -> 3388 rows × 5 columns
