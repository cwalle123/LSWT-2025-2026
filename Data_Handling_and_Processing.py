"""This file is used to import, export, and save all revelant data for easy access in other files"""

####################################################################################################

# External imports
import csv
import re
import pandas as pd

# Internal imports
import Ambient_Conditions as AC

####################################################################################################
"""Functions"""

def convert_raw_txt_to_csv(input_txt_path, output_csv_path):
    """
    Convert raw_Group5_2d.txt or raw_Group5_3d.txt to a CSV file.
    Handles 2-line headers and whitespace-separated columns.
    """

    with open(input_txt_path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]

    # --- Extract and merge the two header rows ---
    header_row1 = re.split(r"\s+", lines[0])
    header_row2 = re.split(r"\s+", lines[1])

    # Merge headers: use row1 name, and if row2 has a unit/value, append in parentheses
    header = []
    for h1, h2 in zip(header_row1, header_row2):
        if h2 in ["/", "-", "_"]:
            header.append(h1)
        else:
            header.append(f"{h1} ({h2})")

    # --- Parse the data rows ---
    data_rows = []
    for line in lines[2:]:
        row = re.split(r"\s+", line)
        data_rows.append(row)

    # --- Write CSV ---
    with open(output_csv_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)
        writer.writerows(data_rows)

    print(f"Converted: {input_txt_path} → {output_csv_path}")

def load_data(dimension: str):
    """
    Load the Group5 dataset as a pandas DataFrame.
    dimension: "2D" or "3D"
    """
    
    # Normalize input (allow lowercase, etc.)
    dim = dimension.strip().upper()

    if dim == "2D":
        csv_path = "Data/Processed_Data/processed_data_Group5_2d.csv"
    elif dim == "3D":
        csv_path = "Data/Processed_Data/processed_data_Group5_3d.csv"
    else:
        raise ValueError('dimension must be "2D" or "3D"')

    df = pd.read_csv(csv_path)
    return df

####################################################################################################
"""Run this file"""

def main():
    two_dimensional_data = load_data("2D")
    angle_of_attack = two_dimensional_data["Alpha (degrees)"]
    print(angle_of_attack)

if __name__ == "__main__":
    main() # makes sure this only runs if you run *this* file, not if this file is imported somewhere else
