## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(foreach)
# library(remeloc)
devtools::load_all()
library(tidyverse)
d <- 2
manifold <- spaceHyperbolic(d, model = 'ball', r = 2)
coef.shrink <- 0.75 # to avoid being too close to "boundary"
local.reach <- 0.2

## ----fig.width=4, fig.height=4------------------------------------------------
set.seed(1)
obsv.coord <- manifold$genPnt(10^5)
idx.edge <- genGraph(10^6, obsv.coord, local.reach = local.reach, exact = F)
plt.dat <- plotGraph(obsv.coord, idx.edge, max.n.edge = 250) +
  coord_fixed(ratio = 1)
plt.dat

## -----------------------------------------------------------------------------
# squared distance
dat <- tibble(sqdist = manifold$dist(
  obsv.coord[idx.edge[, 1], ], obsv.coord[idx.edge[, 2], ]
) ^ 2)
dat <- dat %>% mutate(
  noise = rnorm(length(sqdist), 0, sqrt(sqdist * 0.1)),
  y = pmax(0, sqdist + noise)
) %>% bind_cols(idx.edge %>% as.data.frame)
names(dat) <- c('sqdist', 'noise', 'y', sprintf('p%s', seq(2)))
dat %>% head

## -----------------------------------------------------------------------------
fit <- fitMetric(
  y ~ p1 : p2, 
  data = dat, coord = obsv.coord, 
  optns = list(n.local = 10^4)
)

## -----------------------------------------------------------------------------
target <- manifold$genGrid(128) * coef.shrink
# parallel estimation
system.time({
  N.CORES <- 20
  doParallel::registerDoParallel(cores = N.CORES)
  st.time <- proc.time()
  esti.metric <- foreach(
    w.target = split.data.frame(
      target, # %>% head(2 * N.CORES), 
      cut(seq(nrow(target)), N.CORES)
    ),
    .export = 'fit', .packages = 'remeloc', .errorhandling = 'pass'
  ) %dopar% {
    estiMetric(w.target, fit)
  }
  esti.metric <- unlist(esti.metric, recursive = F)
})

## ---- fig.width=6.5, fig.height=6.5-------------------------------------------
df.esti.metric <- getMetricDF(target, esti.metric)
df.true.metric <- getMetricDF(target, manifold$metric)
map.esti <- df.esti.metric %>% mutate(what = 'esti') %>% 
  bind_rows(df.true.metric %>% mutate(what = 'true')) %>% 
  bind_rows(
    left_join(df.true.metric %>% rename('true' = 'value'), df.esti.metric) %>% 
      mutate(value = true - value, what = 'err') %>% select(-all_of('true'))
  ) %>% 
  filter(component != 'dx1dx2', what != 'err') %>%
  # filter(what == 'err') %>%
  ggplot() +
  geom_raster(aes(x = x1, y = x2, fill = value)) +
  scale_fill_gradientn(
    colors = hcl.colors(10, 'Zissou 1'),
    trans = scales::pseudo_log_trans(base = 10)
  ) +
  facet_wrap(what ~ component, labeller = label_both) +
  coord_fixed(ratio = 1)
map.esti

## ----fig.width=6.5, fig.height=4----------------------------------------------
approx.metric <- approxfunMetric(target, esti.metric)
geo.f.true <- getGeodesic(manifold$metric, d = d)
geo.f.esti <- getGeodesic(function(x){estiMetric(x, fit)[[1]]}, d = d)
geo.f.approx <- getGeodesic(approx.metric, d = d)
df.geo <- bind_rows(
  geo.f.true(rep(0.75, d), c(-1, 0)) %>% as.data.frame %>% as_tibble %>%
    mutate(what = 'true'),
  geo.f.esti(rep(0.75, d), c(-1, 0)) %>% as.data.frame %>% as_tibble %>% 
    mutate(what = 'estimated'),
  geo.f.approx(rep(0.75, d), c(-1, 0)) %>% as.data.frame %>% as_tibble %>% 
    mutate(what = 'approx')
)
df.geo %>% ggplot() +
  geom_path(aes(x = x1, y = x2, color = what)) +
  coord_fixed()

