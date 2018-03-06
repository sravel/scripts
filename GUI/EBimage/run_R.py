from rpy2.robjects.packages import importr
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage


# import R's "base" package
base = importr('base')

with open("fonctions_apprentissage.r", "r", encoding="utf-8") as apprentissageRopen:
	apprentissage = "".join(apprentissageRopen.readlines())
print(apprentissage)


apprentissage = SignatureTranslatedAnonymousPackage(apprentissage, "apprentissage")


path_sample = "/media/sebastien/Bayer/ScriptsSEB/scripts/GUI/EBimage/AnalyseImagesV2/Samples/5583"

print(dir(apprentissage))

print(apprentissage.apprentissage(path_sample))
