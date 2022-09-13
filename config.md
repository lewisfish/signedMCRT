# Config file settings

The below sections describe the tables (dictionaries) that are able to be defined for SignedMCRT

## source

This table defines the parameters for the light source used in the simulation it can have the following:

| Parameter | Type | Options | Default | Notes |
|:---------:|:----:|:-------:|:-------:|:----:|
| name | string | point, circular, uniform, pencil, annulus, focus | point | - |
| nphotons | integer | - | 1000000 | - |
| position | float array size 3 | - | [0.0, 0.0, 0.0] | Default value only set for point source type|
| direction | float array size 3 or string | - | -z | String type applies to all source types bar: Uniform and circular |
| point1 | float array size 3 | - | [-1.0, -1.0, -1.0] | Used by uniform source only to set location and size of source |
| point2 | float array size 3 | - | [2.0, 0.0, 0.0] | See Above |
| point3 | float array size 3 | - | [0.0, 2.0, 0.0] | See Above |
| Radius | float | - | 0.5 | Only used by circular source |

## grid

| Parameter | Type | Options |
|:---------:|:----:|:-------:|
| nxg | integer | - |
| nyg | integer | - |
| nyg | integer | - |
| xmax | float | - |
| ymax | float | - |
| zmax | float | - |
| units | string | - |

## geometry

| Parameter | Type | Options |
|:---------:|:----:|:-------:|
| geom_name | string | - |
| tau | float | - |
| num_spheres | integer | - |
| musb | float | - |
| muab | float | - |
| musc | float | - |
| muac | float | - |
| hgg | float | - |


## detectors

| Parameter | Type | Options |
|:---------:|:----:|:-------:|
| type | string | annulus, circle |
| position | float array size 3 | - |
| radius1 | float | - |
| radius2 | float | - |
| layer | integer | - |
| nbins | integer | - |
| maxval | float | - |

## output

| Parameter | Type | Options |
|:---------:|:----:|:-------:|
| fluence | string | - |
| render | string | - |
| render_geom | boolean | - |
| render_size | integer array size 3 | - |
| overwrite | boolean | - |

## simulation

| Parameter | Type | Options |
|:---------:|:----:|:-------:|
| iseed | integer | - |
| tev | boolean | - |
