import pandas as pd

# Criar função para filtrar o ano e a cidade desejada
def filtrar_dados(data, city, year):
    # Criando uma coluna no df que contenha o ano  
    data["year"] = pd.to_datetime(data["datahora"]).dt.year
    #Filtrando pela cidade desejada
    dt_filtrada_city = data[data["nome_munic"] == city].copy()
    #Pegando a quantidade de casos de covid acumulada no último dia do ano anterior
        # Se o ano for maior que 2020 (já que em 2019 não teria casos)
    if year>2020:
        filtro_casos = dt_filtrada_city["datahora"] == f"{year-1}-12-31"
        numb_correct = (dt_filtrada_city[filtro_casos]["casos"]
                        .iloc[-1])
    else:
        numb_correct = 0
    # Filtrando o ano desejado
    dt_filtrada = dt_filtrada_city[dt_filtrada_city["year"] == year].copy()
    # Realizando a correção, caso seja necessária
    dt_filtrada["casos_sub"] = dt_filtrada["casos"] - numb_correct
    # Filtrando zeros 
    dt_filtrada = dt_filtrada[dt_filtrada["casos_sub"] != 0].copy()
    # Realizando a média dos casos acumulados dos últimos 7 dias
    dt_filtrada["casos_sub"] = (dt_filtrada["casos_sub"]
                                .rolling(window=7)
                                .mean()
                                .bfill())
    return dt_filtrada["casos_sub"]

data = pd.read_csv("../data/dados_covid_sp.csv", sep = ";")

cidades = ["São Paulo", "Bauru"]
anos = [2020, 2021, 2022]

for i in cidades:
    for k in anos:
        (filtrar_dados(data, i, k)
         .to_csv(f"../data/filtered_data/{i}_casos_mm7d_{k}", 
                 index= False))