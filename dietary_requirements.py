import pandas as pd
import enum
import numpy as np

class ActivityLevel(enum.IntEnum):
    SEDENTARY = 0
    MODERATE = 1
    ACTIVE = 2
    VERY_ACTIVE = 3

ActivityLevel_docs = {
    ActivityLevel.SEDENTARY: "light exercise",
    ActivityLevel.MODERATE: "at least 30min of moderate or vigorous exercise every day",
    ActivityLevel.ACTIVE: "at least 1h of vigorous exercise every day",
    ActivityLevel.VERY_ACTIVE: "several hours of vigorous exercise every day",
}


class Sex(enum.Enum):
    FEMALE = 1
    MALE = 2


class AgeGroup(enum.Enum):
    TODDLER = 1
    CHILD = 2
    ADULT = 3


class Nutrient(enum.IntEnum):
    Energy = 1008
    Protein = 1003
    Fat = 1004
    Carb = 1005
    Fiber = 1079
    Linoleic_Acid = 1269
    A_Linoleic_Acid = 1270
    Alcohol = 1018
    Trans_fats = 1257
    Saturated_fats = 1258
    Cholesterol = 1253
    Calcium = 1087
    Chromium = 1096
    Copper = 1098
    Fluoride = 1099
    Iodine = 1100
    Iron = 1089
    Magnesium = 1090
    Manganese = 1101
    Molybdenum = 1102
    Phosphorus = 1091
    Selenium = 1103
    Zinc = 1095
    Potassium = 1092
    Sodium = 1093
    Chloride = 1088  # Chlorine in nutrient table
    Boron = 1137
    Nickel = 1146
    Vitamin_A = 1106
    Vitamin_C = 1162
    Vitamin_D = 1114
    Vitamin_E = 1158
    Vitamin_K = 1185  # only phylloquinone, USDA seems to use that
    Thiamin = 1165
    Riboflavin = 1166
    Niacin = 1167
    Vitamin_B6 = 1175
    Vitamin_B12 = 1178
    Folate = 1177
    Pantothenic_Acid = 1170
    Biotin = 1176
    Choline = 1180


    @classmethod
    def reverse(klass):
        return dict([(int(v), k) for k, v in klass.__members__.items()])


# Table 2 NAS report
ime_physical_activity = {
    (AgeGroup.ADULT, Sex.FEMALE): [1, 1.12, 1.27, 1.45],
    (AgeGroup.ADULT, Sex.MALE): [1, 1.11, 1.25, 1.48],
    (AgeGroup.CHILD, Sex.FEMALE): [1, 1.16, 1.31, 1.56],
    (AgeGroup.CHILD, Sex.MALE): [1, 1.13, 1.26, 1.42],
    (AgeGroup.TODDLER, Sex.FEMALE): [1, 1, 1, 1],
    (AgeGroup.TODDLER, Sex.MALE): [1, 1, 1, 1],  # Same as female
}


# Table 1 NAS report
ime_coeffs = {
    (AgeGroup.ADULT, Sex.FEMALE): [354, 6.91, 9.36, 726],
    (AgeGroup.ADULT, Sex.MALE): [662, 9.53, 15.91, 539.6],
    (AgeGroup.CHILD, Sex.FEMALE): [135.3, 30.8, 10, 934],
    (AgeGroup.CHILD, Sex.MALE): [88.5, 61.9, 26.7, 903],
    (AgeGroup.TODDLER, Sex.FEMALE): [-80, 0, 89, 0],
    (AgeGroup.TODDLER, Sex.MALE): [-80, 0, 89, 0],  # same as female
}

MONTH = 1/12
# Account for energy for "tissue deposition", which I guess is growth. It
# varies by age and is the same across sexes.
# Table 1 NAS report
ime_deposition = [
    (3*MONTH, 175),
    (6*MONTH, 56),
    (1, 22),
    (8, 20),
    (18, 25),
    (np.inf, 0),
]

def cutoff_decision(branches, f):
    for cutoff, decision in branches:
        if f(cutoff):
            return decision
    raise ValueError(f"No branch for value {value} in {branches}")

def eer_kcal(age, sex, weight, height, activity_level):
    age_group = cutoff_decision([
        (2, AgeGroup.TODDLER),
        (18, AgeGroup.CHILD),
        (np.inf, AgeGroup.ADULT)
    ], lambda a: age <= a)
    key = (age_group, sex)

    c = ime_coeffs[key]
    PA = ime_physical_activity[key][activity_level]

    eer = c[0] - c[1]*age + PA*(c[2]*weight + c[3]*height)
    deposition = cutoff_decision(ime_deposition, lambda a: age <= a)
    return float(eer + deposition)


def pri_protein(age, sex, weight):
    df = pd.read_csv("eu-tables/protein-dri.csv", index_col="Age")
    df = df.sort_index(ascending=False)
    key = ('Females' if sex == Sex.FEMALE else 'Males')
    g_prot_kg_subject = cutoff_decision(
        df[key].items(),
        lambda a: age >= a)
    return g_prot_kg_subject * weight


def calc_amdr_percent(age):
    "Acceptable Macronutrient Distribution Range"
    fat_range, protein_range = cutoff_decision([
        (4, ((30, 40), (5, 20))),
        (19, ((25, 35), (10, 30))),
        (np.inf, ((20, 35), (10, 35))),
    ], lambda a: age < a)
    return {Nutrient.Fat: fat_range,
            Nutrient.Linoleic_Acid: (5, 10),
            Nutrient.A_Linoleic_Acid: (0.6, 1.2),
            Nutrient.Carb: (45, 65),    # Carbohydrate
            Nutrient.Protein: protein_range}


# Table 3 NAS report
kcal_g_macro = {
    Nutrient.Carb: 4.,
    Nutrient.Fat: 9.,
    Nutrient.Protein: 4.,
    Nutrient.Alcohol: 7.,
}


def calc_amdr_g(amdr_percent, energy_kcal):
    d = {}
    for k, (percent_lb, percent_ub) in amdr_percent.items():
        k_macro = k
        if k in [Nutrient.Linoleic_Acid, Nutrient.A_Linoleic_Acid]:
            k_macro = Nutrient.Fat

        macro_g_lb, macro_g_ub = map(
            lambda v: (energy_kcal * v)/(100 * kcal_g_macro[k_macro]),
            (percent_lb, percent_ub))
        d[k] = (macro_g_lb, macro_g_ub)
    return d


def convert_units(df, units_from, units_to):
    init_len = len(df.dropna())
    df = pd.concat([df, units_from, units_to], axis=1).dropna()

    unit_div = {
        'G': 1,
        'MG': 1000,
        'MG_ATE': 1000,  # MG α-tocopherol equivalent, vit E
        'UG': 1000000,
    }
    def _convert(_x):
        v, u_from, u_to = _x
        return float(v) * unit_div[u_to] / unit_div[u_from]

    df = df.apply(_convert, axis=1, raw=True)
    assert len(df) == init_len, "We lost rows"
    return df


def make_nuts_indices(df):
    init_len = len(df)
    a = pd.Series(Nutrient.reverse(), name='names')
    df = (pd.merge(a, df, left_on='names', right_index=True)
          .drop(['names'], axis=1))
    assert len(df) == init_len, "We lost rows"
    return df


def read_elem(age, sex, df, name, nutrient):
    sex_str = 'Female' if sex == Sex.FEMALE else 'Male'
    df = df.fillna({'Type': sex_str})

    by_sex = df[df.Type == sex_str]
    i = cutoff_decision(
        list(zip(
            map(float, by_sex.Age),
            by_sex.index)),
        lambda a: age <= a)
    reqs = by_sex.loc[i].drop(['Type', 'Age'])
    reqs = make_nuts_indices(reqs)

    units_from = df[df.Type == 'Units'].transpose().drop(['Type', 'Age'])
    units_from = make_nuts_indices(units_from)

    units_to = nutrient.set_index('id').unit_name

    df = convert_units(reqs, units_from, units_to)
    df.name = name
    return df

def elements(age, sex, nutrient):
    dri = read_elem(age, sex, pd.read_csv("us-tables/elem-dri.csv"),
                    "lower", nutrient)
    ul = pd.read_csv("us-tables/elem-ul.csv").drop(
        ['Vanadium', 'Sulfate', 'Silicon'], axis=1)
    ul = read_elem(age, sex, ul, "upper", nutrient)
    return pd.concat([dri, ul], axis=1)


def vitamins(age, sex, nutrient):
    dri = read_elem(age, sex, pd.read_csv("us-tables/vit-dri.csv"), "lower",
                    nutrient)
    ul = read_elem(age, sex, pd.read_csv("us-tables/vit-ul.csv"), "upper",
                   nutrient)
    return pd.concat([dri, ul], axis=1)


def dietary_requirements(age, sex, weight_kg, height_m, activity_level):
    energy_kcal = eer_kcal(age, sex, weight_kg, height_m, activity_level)
    amdr_percent = calc_amdr_percent(age)
    amdr_g = calc_amdr_g(amdr_percent, energy_kcal)
    df = pd.DataFrame.from_dict(amdr_g, orient='index', columns=('lower', 'upper'))

    df.loc[Nutrient.Fiber, "lower"] = 14 * energy_kcal / 1000

    df.loc[Nutrient.Energy, "lower"] = energy_kcal * 0.98
    df.loc[Nutrient.Energy, "upper"] = energy_kcal * 1.02

    # Cholesterol as low as possible
    df.loc[Nutrient.Cholesterol, "cost"] = 0.01
    # Trans fat, saturated fat ALAP
    df.loc[Nutrient.Trans_fats, "cost"] = 0.2
    df.loc[Nutrient.Saturated_fats, "cost"] = 0.1
    # Added sugars ≤ 25%

    # Amend protein lower bound if necessary
    idx = (Nutrient.Protein, "lower")
    df.loc[idx] = max(df.loc[idx], pri_protein(age, sex, weight_kg))

    nutrient = pd.read_csv("FoodData_Central_csv_2019-12-17/nutrient.csv.gz")
    elem = elements(age, sex, nutrient)
    vit = vitamins(age, sex, nutrient)
    return pd.concat([df, elem, vit], axis=0, sort=True)



#if __name__ == '__main__':
df = dietary_requirements(24, Sex.MALE, 88.8, 1.98, ActivityLevel.MODERATE)

nutrient = pd.read_csv("FoodData_Central_csv_2019-12-17/nutrient.csv.gz")

def check_Nutrients():
    a = pd.Series(Nutrient.reverse(), name='bla')
    table = pd.merge(nutrient, a, left_on='id', right_index=True)[['id', 'bla', 'name']]
    print(table)


with_nutrient = pd.merge(nutrient[['id', 'name']], df, left_on='id', right_index=True)
a = pd.Series(Nutrient.reverse(), name='names')

print(pd.concat([a, with_nutrient.set_index('id')], axis=1))
