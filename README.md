

# Cosmological analyses made easy

### Installation

Installation should be quite easy. Either download the zip or clone the repo, and inside the folder do

`python3 pip install -e .`

If this doesn't work, please send me the error message that you're getting!

### Basic Usage

The usage should in principle be relatively straightforward. The main functionality is `load` to load a collection (or single) folder

    import liquidcosmo as lc

    fo = lc.load("foldername1","foldername2",...)

    spp = fo[parameters].plot_getdist()

    spp.export("file.pdf")

## Documentation

(There will be future updates to this section, for now it's listing just the most important functions)

	load(path, [burnin_threshold, timeout, tag])

Here `path` is the path to the chain folder containing the chain files, or the path compared to a `chains` folder in the same directory. You can also give individual file paths.

`burnin_threshold` is the MontePython way of removing burnin. Here the first items of the chain are all removed until a loglike is reached that is at least as small as the bestfit loglike with a offset of `+burnin_threshold`. (the default is `+3` for MontePython consistency)

`timeout` allows you to wait for a longer for chains to be loaded from your file system. The default is `60` seconds.

Finally, `tag` is a possibility to give the folder a name that will be put in the legend of a possible plot. If none is given, then it takes the foldername instead. It can always be overwritten by using the `legend_labels`option of `plot_getdist`, see below.

The result of this function will be a `foldercollection` if you loaded multiple folders, or a `folder`if you loaded a single folder. The two act mostly the same for most purposes. See below for the documentation.

### Handling `folder`/`foldercollection` objects

#### Modifying the object

Let's call the result of the `load` function simply as `fo`. We have the following things that we can do

	fo[x]

Here x can be a large number of things, it can be:
 * A comma-seperated list of parameter names like `fo['name1','name2']` or a python list or tuple of the same `fo[list]` with `list=['name1','name2']`.
 * A single integer number `i`. For a `foldercollection` this will give you the `i-th` folder, while for a `folder` it will give you the `i-th` element of the chain.
 *  For a `foldercollection`, you can also pass a python slice, like `a:b` to return folders from `a` to `b`.
 * You can also give a list of booleans (a mask), which will allow you to select only the items in the folder which suffice the boolean condition (have `True` at the correct locations). For example, it is possible to give an array like `fo[mask]` with `mask = fo['omega_b']>0.022` (see below) to select only the elements for which the `omega_b` parameter is bigger than `0.022`.

There are a bunch of other similar and useful functions as well:

	fo.cut([first, last, thin])

You can cut the folder items to only keep part of the chain. `first` is a percentage (number in `[0..1]`) that is the percentage of chain that is removed from the start of the chain. Similarly, `last` is also a percentage, but this time specifying the upper bound (until which percentage point are kept). Finally `thin` is an integer number, where only ever `thin-th` point is kept.

	fo.set_range(parname, [lower,upper,destructive])

The `parname` specifies the parameter name, the `lower`is the lower range that should be applied to the parameter, while `upper` is the upper range applied to the parameter. Finally, if `destructive` is not used, the folder is not actually reduced (the points outside the range are kept, despite being outside of the new range).

	fo.to_individual()

Splits the `folder` object into a `foldercollection`, containing the separate chains in the folder.

	fo.merge([basechain])

The merge function merges a `foldercollection` into a single `folder` object (undoing `to_individual()`, but it can also merge independent chains)

#### Getting information from the objects

There are also a bunch of functions to retrieve information from the `folder`/`foldercollection` objects.

	fo.bestfit()

Gives the bestfit or bestfits.

	fo.names

Gives the names of the parameters in the chain. To get just the parameters, you can use `fo.names[2:]`.

	fo.mean([parnames,asdict])
	fo.cov([parnames])
	fo.std([asdict])

This gives the mean, covariance matrix, or standard deviation. The `parnames` is a list of parameter names (or a string of a single parameter). The `asdict` parameter allows you to retrieve the output as a dictionary with the parameter names and then the corresponding values.

	fo.set_texname(parname, texname)

Used to set the tex name for a parameter given by `parname` as `texname` -- The latter should be render-able with latex.

	fo.get_bounds()

This gives the ranges of the parameters as a dictionary. Similarly `get_range(par)` gives the range for a single parameter `par`.

	fo.write(fname, [codetype])

The `fname` is writing the current folder object (after modifications) in the same format as the code `codetype`.

Getting `samples` returns an array of the samples, `logfile`the arguments of the file containing the properties of the chain, similarly `cosmoargs` gives only the fixed arguments.

	fo.to_class_dict()
	fo.to_class_ini()

Can be used to turn the chain object into a python dictionary, or a `.ini` file to be used by `class`. The second is only possible if there's only one point in the chain left, such as using the `fo.bestfit.to_class_ini()`.

	fo.constraint([parnames])

This gives the constraints in a special notation, for each parameter in `parnames` (either string or list of strings). It gives `[[lower, upper], type]`, where `lower` is the lower bound (or lower range), `upper` is the upper bound (or upper range), and `type` can be `unconstrained` if there is no constraint, `>` if there is only a lower bound, and `<` if there is only an upper bound, and `+-` if there is a lower and upper bound (i.e. sigma errors).

	fo.texconstraint([parnames,withdollar,withname])

This gives the constraints for the parameters specified in `paramnames` (single string or list thereof) well formatted for inserting e.g. into a latex document.  `withdollar` says whether a `$` prefix/suffix should be included, and `withname` says whether the name should be included in the resulting string.

	fo.gelman([parnames,subdivisions])

Gives the Gelman-Rubin criterion for the parameters specified by `paramnames` (single list or list thereof). `subdivisions` can be used to additionally subdivide the array.

	fo.max_gelman([subdivisions])

Gives the maximum Gelman-Rubin criterion for all possible parameter directions (which are specified by the eigenvalues of the covariance matrix). Can be subdivided with `subdivisions`, see `gelman` above.

#### Plotting

There is also a functions convert the `folder` or `foldercollection` objects using the [`getdist`](https://getdist.readthedocs.io/en/latest/index.html)  tool, and to plot it.

	fo.to_getdist()

Converts the folder into a `MCSamples` object from `getdist`

	fo.plot_getdist([ax, colors, alphas, add_point, show, contours, **kwargs])

This plots the object by first converting them to an `MCSamples` object , and the invoking the plot option. The `alphas` specifies the transparency as an `alpha` values for each contour. The `add_point` adds a new point (as vertical/horizontal lines) to the plot. `show` invokes `matplotlib.pyplot.show` at the end of the function. The `contours` says how many contours should be drawn, either as an integer (number of sigma intervals) or as a list of percentages. Finally, `kwargs` is a list of additional keyword arguments that will generally be passed to the `triangle_plot` or `rectangle_plot` function of `getdist`, except for a few exceptions
 * `legend_labels` are the labels that should be in the legend of the plot.
 * `analysis_settings` are `getdist` analysis settings that are passed to the `updateSettings` function of `getdist`, see [this link here](https://getdist.readthedocs.io/en/latest/analysis_settings.html).
 * `linestyle` is a list of linestyles (or a single one) that should be applied to the lines (in the 1D posteriors)
 * Similarly, `line_args` is a list of dictionaries of the arguments for each line for each folder that is plotted.
 * If the option `rectangle` is included, it should be a dictionary containing `x` and `y` for the parameters to draw on the x-axis/y-axis, each one a list of parameter names, for example:
	 `fo.plot_getdist(rectangle={'x':['a','b'],'y':['c','d','e']})`
