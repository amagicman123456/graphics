<div align="center">

## the project (made by ryu co)

this is a raytracer, currently able to render spheres, polygons, and tori

![image](/images/line-torus.PNG)

but idk what im doing

btw outdated/yes_cursed.exe is way too cursed

![image](/images/cursed.PNG)

</div>

build with

```
x86_64-w64-mingw32-g++ yes.cpp [options] -static -lgdi32 -fconcepts-ts -Ofast -o yes
```

options:

* -DNO_DOUGHNUT: no doughnut in sight :(

* -Dtri_specialization: use an algorithm specifically for triangles, idk if the performance changes

* -Dfirst_algorithm: use the first of the triangle algorithms

* -DQUAD: render a quadrilateral instead of a triangle, i've neglected this so idk
