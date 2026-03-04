# Importe the Pandas library for filter and manipulate COVID data from the municipality of São Paulo
import pandas as pd

def data_filter(data:pd.DataFrame, city:str, year:int):

    # Create a new column "year" on a dataframe extracting the year of the column "datahora"    
    data["year"] = pd.to_datetime(data["datahora"]).dt.year
    
    # Filter the selected city and create a copy of the DataFrame
    filtered_city_data = data[data["nome_munic"] == city].copy()

    # Adjust cumulative baseline using the previous year's final value (if applicable)
    if year > 2020:
        # Create a filter for December 31st of the previous year
        previous_year_end_filter = filtered_city_data["datahora"] == f"{year-1}-12-31"
        # Retrieve the accumulated number of cases on December 31st of the previous year 
        numb_correction = (
            filtered_city_data[previous_year_end_filter]["casos"]
            .iloc[-1]
        )
    # If the year is 2020, no correction is needed
    else:
        numb_correction = 0

    # Filter the selected year and create a copy
    filtered_data = filtered_city_data[filtered_city_data["year"] == year].copy()
    
    # Convert cumulative cases to yearly cases using a correction 
    filtered_data["cases_7-dma"] = filtered_data["casos"] - numb_correction

    # Remove initial zeros from the dataset 
    filtered_data = filtered_data[filtered_data["cases_7-dma"] != 0].copy()
    
    # Compute 7-day moving average of the cumulative cases and Backfill initial NaN values
    filtered_data["cases_7-dma"] = (
        filtered_data["cases_7-dma"]
        .rolling(window=7)
        .mean()
        .bfill()
    )

    # Return 7-day moving average for the year and municipality selected
    return filtered_data["cases_7-dma"]

# Load COVID-19 case data for the state of São Paulo to begin the data filtering
data = pd.read_csv("../data/dados_covid_sp.csv", sep=";")

# Define municipalities to filter COVID-19 data for São Paulo and Bauru (2020 - 2022) 
municipalities = ["São Paulo", "Bauru"]
years = [2020, 2021, 2022]

# Iterate over municipalities and years to export COVID-19 data (2020-2022) 
for municipality in municipalities:
    for year in years:
        filtered_df = data_filter(data, municipality, year)
        filtered_df.to_csv(
            f"../data/filtered_data/{municipality}_casos_mm7d_{year}", 
            index=False
        )
        