# %%
# Importe the Pandas library for filter and manipulate COVID data from the municipality of São Paulo
import pandas as pd

# Create a function for date and year filtering 
def filtrar_dados(data:pd.DataFrame, city:str, year:int):

    # Create a new column "year" on a dataframe extracting the year of the column "datahora"    
    data["year"] = pd.to_datetime(data["datahora"]).dt.year
    
    # Filter the selected city and create a copy of the DataFrame
    dt_filtrada_city = data[data["nome_munic"] == city].copy()

    # Previous year based on the selected year
    # Use an if/else statement for saving the total number of accumulated cases, 
    # on December 31st of the previous year 
    # This only applies if the year is greater than 2020,
    # because on 2019 we didn't have any cases of COVID-19
    if year > 2020:
        # Create a filter for December 31st of the previous year
        filtro_casos = dt_filtrada_city["datahora"] == f"{year-1}-12-31"
        # Retrieve the accumulated number of cases on December 31st of the previous year 
        numb_correction = (dt_filtrada_city[filtro_casos]["casos"]
                        .iloc[-1])
    # If the year is 2020, no correction is needed
    else:
        numb_correction = 0

    # Filter the selected year and create a copy
    dt_filtrada = dt_filtrada_city[dt_filtrada_city["year"] == year].copy()
    
    # Adjust cumulative cases 
    dt_filtrada["cases_7-dma"] = dt_filtrada["casos"] - numb_correction

    
    # Filter zeros and create a copy 
    dt_filtrada = dt_filtrada[dt_filtrada["cases_7-dma"] != 0].copy()
    
    # Calculate the 7-day moving average of the cumulative cases and
    # Backfill missing values using the next valid observation 
    # to avoid NaNs at the beginning of the series   
    dt_filtrada["cases_7-dma"] = (
        dt_filtrada["cases_7-dma"]
        .rolling(window=7)
        .mean()
        .bfill()
    )

    # Return a Series with the moving avarage 
    return dt_filtrada["cases_7-dma"]

data = pd.read_csv("../data/dados_covid_sp.csv", sep = ";")

cidades = ["São Paulo", "Bauru"]
anos = [2020, 2021, 2022]

for i in cidades:
    for k in anos:
        (filtrar_dados(data, i, k)
         .to_csv(f"../data/filtered_data/{i}_casos_mm7d_{k}", 
                 index= False))
        