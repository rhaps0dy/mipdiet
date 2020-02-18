import pandas as pd
import os
import gzip

path_to_db = "./FoodData_Central_csv_2019-12-17/"

nutrient = pd.read_csv(os.path.join(path_to_db, "nutrient.csv"))
food_nutrient = pd.read_csv(os.path.join(path_to_db, "food_nutrient.csv"))
food = pd.read_csv(os.path.join(path_to_db, "food.csv"))

# This is the future table of foods with nutrients
food = food[["fdc_id", "description"]].set_index("fdc_id")
assert not food.index.duplicated().any()
nutrient = nutrient.set_index("id")
assert not nutrient.index.duplicated().any()

# make nutrients the columns, and foods the rows
fused_nutrients = food_nutrient.pivot(
    index="fdc_id", columns="nutrient_id", values="amount")
# "data_points", "min", "max", "median"])

df = food.merge(fused_nutrients, how="inner", left_index=True, right_index=True)
df.to_csv("foods.csv.gz", compression="gzip")

s = "\n".join(
    f"(\"{desc}\" . \"{fdc_id}\")"
    for fdc_id, desc in df['description'].iteritems())
s = f"(quote ({s}))"

with gzip.open("food-names.el.gz", "wt", encoding="UTF-8") as f:
    f.write(s)

