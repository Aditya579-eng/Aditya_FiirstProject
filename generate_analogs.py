import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import pandas as pd
import random

# Function to calculate molecular properties
def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        mw = round(Descriptors.MolWt(mol), 2)  # Molecular weight
        logp = round(Descriptors.MolLogP(mol), 2)  # Lipophilicity
        tpsa = round(rdMolDescriptors.CalcTPSA(mol), 2)  # Polar Surface Area
        hbd = rdMolDescriptors.CalcNumHBD(mol)  # H-bond donors
        hba = rdMolDescriptors.CalcNumHBA(mol)  # H-bond acceptors
        return mw, logp, tpsa, hbd, hba
    return None, None, None, None, None

# Function to generate analogs
def generate_analogs(smiles, num_analogs=10):
    mol = Chem.MolFromSmiles(smiles)
    analogs = []

    if mol:
        for _ in range(num_analogs):
            new_mol = Chem.RWMol(mol)

            # Randomly mutate the structure by adding atoms (C, N, O)
            atom_idx = random.choice(range(new_mol.GetNumAtoms()))
            new_atom = random.choice([6, 7, 8])  # Carbon, Nitrogen, Oxygen
            new_mol.AddAtom(Chem.Atom(new_atom))
            new_mol.AddBond(atom_idx, new_mol.GetNumAtoms() - 1, Chem.BondType.SINGLE)

            # Convert back to SMILES
            new_smiles = Chem.MolToSmiles(new_mol)
            mw, logp, tpsa, hbd, hba = calculate_properties(new_smiles)

            analogs.append([new_smiles, mw, logp, tpsa, hbd, hba])

    return analogs if analogs else [["Invalid SMILES", None, None, None, None, None]]

# Streamlit Web App UI
st.set_page_config(page_title="Phytochemical Analog Generator", page_icon="🧪", layout="wide")
st.title("🧪 Phytochemical Analog Generator")
st.markdown("Enter a **SMILES structure** to generate modified analogs with molecular properties.")

# **Single SMILES input**
smiles_input = st.text_input("🔹 Enter a SMILES string:")
if st.button("Generate Analogs"):
    if smiles_input:
        analogs = generate_analogs(smiles_input)
        df_single = pd.DataFrame(analogs, columns=["Generated SMILES", "MW", "LogP", "TPSA", "HBD", "HBA"])
        st.write(df_single)
        st.download_button("📥 Download Results", df_single.to_csv(index=False), file_name="analogs.csv")
    else:
        st.error("⚠️ Please enter a valid SMILES string.")

# **Batch Processing: Upload Excel File**
uploaded_file = st.file_uploader("📂 Upload an Excel file with SMILES (Column name: 'SMILES')", type=["xlsx"])
if uploaded_file:
    try:
        df = pd.read_excel(uploaded_file)
        if "SMILES" not in df.columns:
            st.error("⚠️ The uploaded file must contain a column named 'SMILES'.")
        else:
            smiles_list = df["SMILES"].dropna().tolist()
            
            results = []
            for smiles in smiles_list:
                analogs = generate_analogs(smiles)
                for analog in analogs:
                    results.append([smiles] + analog)  # Store original & generated analogs

            df_results = pd.DataFrame(results, columns=["Original SMILES", "Generated SMILES", "MW", "LogP", "TPSA", "HBD", "HBA"])
            st.write(df_results)
            st.download_button("📥 Download Results", df_results.to_csv(index=False), file_name="batch_analogs.csv")
    except Exception as e:
        st.error(f"⚠️ Error processing file: {str(e)}")
