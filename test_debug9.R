# the original with(deriv_out, .find_valleys(x, ...)) in cytokine_cutpoint.R looks like a bug 
# since it uses `deriv_out$x` as the input to .find_valleys, which will just compute the density of uniformly spaced points!

# Let's check cytoUtils again. Maybe `deriv_out` in cytoUtils is NOT just a list of x and y.
# If .deriv_density returns a density object, maybe .find_valleys takes a density object?
