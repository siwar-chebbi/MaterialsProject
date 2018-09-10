from pymatgen import MPRester


api = MPRester("fB610TDF3LSwxiN9")

#entries = api.get_entries({"nelements": {'$lte': 6, '$gte': 1},"elasticity": {'$ne': None}}, property_data=['pretty_formula','elasticity', 'elements'])
covalent = ['B', 'C', 'Si']
ionique = ['N', 'O', 'F', 'P', 'S', 'Cl', 'Sr', 'Br', 'I']
#tableau des proprietes
propsTableau = ['material_id','pretty_formula','elasticity.poisson_ratio','elasticity.G_Reuss', 'elasticity.G_Voigt', 'elasticity.G_Voigt_Reuss_Hill', 'elasticity.K_Reuss', 'elasticity.K_Voigt', 'elasticity.K_Voigt_Reuss_Hill']

critere1 = {"nelements": {'$lte': 6},'elements': {'$in': covalent}, "elasticity": {'$ne': None}, "elasticity.poisson_ratio": {'$lt': 0 },'elasticity.G_Reuss' : {'$gte': 0 }, 'elasticity.G_Voigt': {'$gte': 0 }, 'elasticity.G_Voigt_Reuss_Hill': {'$gte': 0 }, 'elasticity.K_Reuss': {'$gte': 0 }, 'elasticity.K_Voigt': {'$gte': 0 }, 'elasticity.K_Voigt_Reuss_Hill': {'$gte': 0 }}
critere2 = {"nelements": {'$lte': 6},'elements': {'$in': ionique}, "elasticity": {'$ne': None}, "elasticity.poisson_ratio": {'$lt': 0 },'elasticity.G_Reuss' : {'$gte': 0 }, 'elasticity.G_Voigt': {'$gte': 0 }, 'elasticity.G_Voigt_Reuss_Hill': {'$gte': 0 }, 'elasticity.K_Reuss': {'$gte': 0 }, 'elasticity.K_Voigt': {'$gte': 0 }, 'elasticity.K_Voigt_Reuss_Hill': {'$gte': 0 }}


materials1 = api.query(criteria=critere1, properties=propsTableau)
materials2 = api.query(criteria=critere2, properties=propsTableau)

def recup(materials):
    texte=""
    for prop in propsTableau:
        texte = texte + str(prop) + "\t"
    print(texte)

    for material in materials:
        texte = ""
        for prop in propsTableau:
            texte = texte + str(material.get(prop)) + "\t"
        print(texte)


print("\n********************ELEMENTS COVALENTS*******************************\n")
recup(materials1)
print("\n********************ELEMENTS IONIQUES*******************************\n")
recup(materials2)