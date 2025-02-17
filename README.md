# Résolveur de Poisson utilisant des méthodes itératives

Ce projet implémente un résolveur de Poisson utilisant plusieurs méthodes itératives, dont **Jacobi**, **Gauss-Seidel** et **Gradient Conjugué**, pour résoudre un système d'équations linéaires provenant de la discrétisation de l'équation de Poisson sur une grille 2D.

## Compilation

Le code est écrit en **C** et peut être compilé avec différents compilateurs. Le `Makefile` fourni permet une compilation automatique avec les options suivantes :

### Options du compilateur :
- `gcc`
- `icc` (Compilateur Intel)
- `icx` (Compilateur oneAPI d'Intel)
- `clang`
- `aocc` (Compilateur AMD)

### Compilation par défaut

Pour compiler avec les paramètres par défaut, exécutez simplement :

```sh
make poisson
```

### Compilation personnalisée

Vous pouvez spécifier le compilateur et les options en modifiant le `Makefile`. Par exemple, pour compiler avec le support de MKL d'Intel, utilisez :

make CC=icc poisson


### Options d'optimisation

Le programme est compilé avec l'option `-O3` pour optimiser les performances.

### Bibliothèques utilisées

Le code est lié avec les bibliothèques suivantes :
- **Intel MKL** (lors de l'utilisation des compilateurs Intel)
- **Bibliothèque mathématique GNU** (`-lm`) pour les autres compilateurs.

### Nettoyage des fichiers générés

Pour nettoyer les fichiers compilés, exécutez :

```sh
make clean
```

## Utilisation

### Exécution du résolveur de Poisson

Une fois compilé, le résolveur peut être exécuté comme suit :

```sh
./poisson <dim> <max_iter> <threshold>
```

Où :
- **dim** : La dimension de la grille (par exemple, 4 pour une grille 4x4, 10 pour une grille 10x10, etc.).
- **max_iter** : Le nombre maximal d'itérations pour l'exécution du résolveur.
- **threshold** : Le seuil de convergence pour l'erreur résiduelle.

Par exemple, pour exécuter le résolveur sur une grille 4x4 avec un maximum de 150 itérations et un seuil de 0.000001 :

```sh
./poisson 4 150 0.000001
```

### Sortie

Le résolveur affichera les informations suivantes :
- Le nombre d'itérations et l'erreur à la fin de chaque méthode itérative.
- Le résultat du résolveur et l'erreur résiduelle, comparée à la solution attendue.

### Méthodes Implémentées

1. **Méthode de Jacobi** : Une méthode itérative pour résoudre des systèmes linéaires, qui met à jour le vecteur solution en résolvant pour chaque élément de manière individuelle en fonction de l'itération précédente.
   
2. **Méthode de Gauss-Seidel** : Similaire à la méthode de Jacobi, mais avec la modification suivante : la solution mise à jour est utilisée immédiatement dans les calculs suivants au sein de la même itération.
   
3. **Méthode du Gradient Conjugué** : Une méthode itérative plus avancée, adaptée aux grands systèmes clairsemés, qui utilise les gradients pour guider la recherche de la solution.

### Représentation de la matrice

La matrice est stockée sous forme **Compressed Sparse Row (CSR)**, une méthode de stockage efficace en mémoire pour les matrices clairsemées. Seules les valeurs non nulles, leurs indices de colonnes et les pointeurs des lignes sont stockés.

### Vérification de l'erreur

Après chaque itération, le programme vérifie si la solution a convergé en calculant l'erreur résiduelle. Si l'erreur est inférieure au seuil spécifié, le résolveur s'arrête, indiquant qu'une solution a été trouvée. L'erreur est calculée comme suit :

residual_error = ||Ax - b|| / ||b||


Où `A` est la matrice, `x` est le vecteur solution et `b` est le vecteur de droite.