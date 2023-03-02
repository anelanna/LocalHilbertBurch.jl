using Documenter
using LocalHilbertBurch

makedocs(
    sitename = "LocalHilbertBurch",
    format = Documenter.HTML(),
    modules = [LocalHilbertBurch]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/lkastner/LocalHilbertBurch.jl.git"
)
