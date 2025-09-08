## ----eval = FALSE-------------------------------------------------------------
# install.packages("bipartite")

## ----results = "hide", message = FALSE, cache = FALSE-------------------------
library(bipartite)
set.seed(123)

## -----------------------------------------------------------------------------
web <- genweb(N1 = 5, N2 = 6)

## -----------------------------------------------------------------------------
str(web)
web
class(web)

## -----------------------------------------------------------------------------
colnames(web) <- paste("Higher", 1:6)
rownames(web) <- paste("Lower", 1:5)
web

## ----eval = FALSE-------------------------------------------------------------
# ?Safariland

## -----------------------------------------------------------------------------
data(Safariland)
dim(Safariland)
str(Safariland)

## ----eval = FALSE-------------------------------------------------------------
# ?plotweb

## ----fig.show="hold"----------------------------------------------------------
plotweb_deprecated(web)
plotweb_deprecated(Safariland)

## -----------------------------------------------------------------------------
plotweb_deprecated(Safariland, text.rot = 90)

## ----first_plot, cache = FALSE------------------------------------------------
plotweb(web)

## ----new_safariland-----------------------------------------------------------
plotweb(Safariland)

## ----safariland_90------------------------------------------------------------
plotweb(Safariland, srt = 90)

## -----------------------------------------------------------------------------
web2 <- matrix(c(50, 50), ncol = 2)
plotweb(web2, plot_axes = TRUE)

## -----------------------------------------------------------------------------
plotweb(web2, spacing = 0.1, plot_axes = TRUE)

## -----------------------------------------------------------------------------
plotweb(web2, spacing = 0.1, text_size = 2)

## ----fig.height=6, fig.width=6------------------------------------------------
plotweb(Safariland, spacing = "auto", srt = 90, text_size = 0.75)

## ----web_horizontal-----------------------------------------------------------
plotweb(web, horizontal = TRUE)

## ----sorting_reverse----------------------------------------------------------
plotweb(web[rev(rownames(web)), ], horizontal = TRUE)

## ----sorting_all--------------------------------------------------------------
plotweb(Safariland, horizontal = TRUE, sorting = "dec")
plotweb(Safariland, horizontal = TRUE, sorting = "inc")
plotweb(Safariland, horizontal = TRUE, sorting = "ca")

## -----------------------------------------------------------------------------
n_lower_species <- nrow(Safariland)
lower_abundances <- sample(0:100, n_lower_species, replace = TRUE)
names(lower_abundances) <- rownames(Safariland)
print(lower_abundances)

## -----------------------------------------------------------------------------
plotweb(Safariland, sorting = "ca", horizontal = TRUE, 
        lower_abundances = lower_abundances)

## -----------------------------------------------------------------------------
plotweb(Safariland, sorting = "ca", horizontal = TRUE,
        add_lower_abundances = lower_abundances)

## -----------------------------------------------------------------------------
n_higher_species <- ncol(Safariland)
higher_abundances <- sample(0:1000, n_higher_species, replace = TRUE)
names(higher_abundances) <- colnames(Safariland)
plotweb(Safariland, sorting = "ca", horizontal = TRUE,
        higher_abundances = higher_abundances, scaling = "relative")
plotweb(Safariland, sorting = "ca", horizontal = TRUE,
        higher_abundances = higher_abundances, scaling = "absolute")

## -----------------------------------------------------------------------------
plotweb(Safariland, sorting = "ca", horizontal = TRUE,
        add_lower_abundances = lower_abundances, scaling = "absolute")

## -----------------------------------------------------------------------------
plotweb(Safariland, srt = 90, curved_links = TRUE)

## -----------------------------------------------------------------------------
plotweb(Safariland, horizontal = TRUE, curved_links = TRUE)

## -----------------------------------------------------------------------------
plotweb(Safariland, sorting = "ca", horizontal = TRUE, curved_links = TRUE,
        higher_italic = TRUE, lower_italic = TRUE)

## -----------------------------------------------------------------------------
higher_labels <- colnames(Safariland)
names(higher_labels) <- higher_labels
species_name_selector <- lengths(strsplit(higher_labels, " ")) == 2
higher_species_names <- higher_labels[species_name_selector]
higher_labels[higher_species_names] <- lapply(higher_species_names,
                                              function(x) bquote(italic(.(x))))

plotweb(Safariland, sorting = "ca", horizontal = TRUE, curved_links = TRUE,
        higher_labels = higher_labels, lower_italic = TRUE)

## -----------------------------------------------------------------------------
plotweb(Safariland, sorting = "ca", horizontal = TRUE, curved_links = TRUE, 
        higher_labels = higher_labels, lower_italic = TRUE,
        higher_color = "orange", lower_color = "darkgreen")


## -----------------------------------------------------------------------------
plotweb(Safariland, sorting = "ca", horizontal = TRUE, curved_links = TRUE,
        higher_labels = higher_labels, lower_italic = TRUE,
        lower_color = rainbow(nrow(Safariland)))

## -----------------------------------------------------------------------------
lower_color <- rep("black", nrow(Safariland))
names(lower_color) <- rownames(Safariland)
lower_color["Alstroemeria aurea"] <- "orange"
plotweb(Safariland, sorting = "ca", horizontal = TRUE, curved_links = TRUE,
        higher_labels = higher_labels, lower_italic = TRUE,
        lower_color = lower_color)


## -----------------------------------------------------------------------------
plotweb(Safariland, sorting = "ca", horizontal = TRUE, curved_links = TRUE,
        higher_labels = higher_labels, lower_italic = TRUE,
        lower_color = lower_color, link_color = "lower")

## -----------------------------------------------------------------------------
plotweb(Safariland, sorting = "ca", horizontal = TRUE, curved_links = TRUE,
        higher_labels = higher_labels, lower_italic = TRUE,
        higher_color = "orange", lower_color = "darkgreen",
        link_color = "brown")

## -----------------------------------------------------------------------------
lower_add_color <- rep("black", nrow(Safariland))
lower_color <- rep("gray50", nrow(Safariland))
names(lower_add_color) <- rownames(Safariland)
names(lower_color) <- rownames(Safariland)
lower_add_color["Alstroemeria aurea"] <- "darkorange"
lower_color["Alstroemeria aurea"] <- "orange"
plotweb(Safariland, sorting = "ca", horizontal = TRUE, 
        add_lower_abundances = lower_abundances, curved_links = TRUE,
        higher_labels = higher_labels, higher_color = "grey50",
        lower_italic = TRUE, lower_color = lower_color,
        link_color = "lower", lower_add_color = lower_add_color)

## -----------------------------------------------------------------------------
lower_add_color <- rep("black", nrow(Safariland))
lower_color <- rep("gray50", nrow(Safariland))
names(lower_add_color) <- rownames(Safariland)
names(lower_color) <- rownames(Safariland)
lower_add_color["Alstroemeria aurea"] <- "darkorange4"
lower_color["Alstroemeria aurea"] <- "orange"
plotweb(Safariland, sorting = "ca", horizontal = TRUE,
        add_lower_abundances = lower_abundances,
        curved_links = TRUE, higher_labels = higher_labels,
        higher_color = "grey50", lower_italic = TRUE,
        lower_color = lower_color, link_color = "lower",
        link_alpha = 1, higher_border = "black",
        link_border = "black", lower_border = "black",
        lower_add_color = lower_add_color)

