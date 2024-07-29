<div align="center">

## the project (made by ryu co)

this is a raytracer in c++ for versions >= c++14, currently able to render spheres, polygons, and tori

![image](/images/line-torus.PNG)

but idk what im doing

btw outdated/yes_cursed.exe is way too cursed

![image](/images/cursed.PNG)

</div>

build with

```
x86_64-w64-mingw32-g++ yes.cpp [options] -static -lgdi32 `# [no longer needed: -fconcepts-ts]` -Ofast -o yes
```

options:

* -DSMOOTHEN: smoothen out the pixels, gives better fps too

* -DSMOOTHEN_AMOUNT=[value]: determines how smoothened the pixels are

* -DNO_DOUGHNUT: no doughnut in sight :(

* -Dtri_specialization: use an algorithm specifically for triangles, idk if the performance changes

* -Dfirst_algorithm: use the first of the triangle algorithms

* -DQUAD: render a quadrilateral instead of a triangle, i've neglected this so idk

* -DSOUND -lwinmm -lstdc++fs: play sounds when clicking shapes :D

* -DIMAGE: render images on the shapes (sphere turns into earth) :D

* -Dhorizontal_fov=[value]: set the horizontal fov to [value] degrees/radians (depending on if -Din_radians is set)

* -Din_radians: specify that the horizontal fov set is in radians
