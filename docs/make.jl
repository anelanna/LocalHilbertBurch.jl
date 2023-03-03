using Documenter
using LocalHilbertBurch

DocMeta.setdocmeta!(LocalHilbertBurch, :DocTestSetup, :(using LocalHilbertBurch, Oscar); recursive=true)

makedocs(
    sitename = "LocalHilbertBurch",
    format = Documenter.HTML(),
    modules = [LocalHilbertBurch]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/anelanna/LocalHilbertBurch.jl.git"
)
