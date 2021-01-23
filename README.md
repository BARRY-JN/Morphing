# Morphing

Projet permettant de générer une animation de morphing entre 2 mesh (inachevé)

Il implémente l'algorithme de cette étude :
Bogdan Cosmin Mocanu - 3D Mesh Morphing - Intitut National des Télécomunications - 2012

Actuellement, la transformation se fait par un appairage point par point, en fonction de leur distance (pour chaque point d'un modèle, on cherche le point de l'autre modèle le plus proche dans l'espace).

Un algorithme de paramétrisation sphérique a été implémenté pour, à terme, utiliser la technique du métamesh.

![morph1](https://user-images.githubusercontent.com/43220602/105612336-a85f0880-5dbb-11eb-8225-adbbaaba5368.png)
![morph2](https://user-images.githubusercontent.com/43220602/105612344-b57bf780-5dbb-11eb-8cbd-15cc76c2d066.png)
