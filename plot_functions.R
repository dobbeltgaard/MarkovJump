library(reshape2)
library(ggplot2)
library(Rcpp)
library(RcppEigen)
sourceCpp("FUNCS_MJP_with_eigen.cpp")
library(Matrix)


#functions to generate figures
get_pars = function(str, fold_number = 1){read.csv(list.files("estimates", full.names = T)[grepl(str, list.files("estimates"))][fold_number])$par} #get_pars("gerlang_relax_FALSE_exp_exp_all_no_warp")

trans_dist_fig <- function(m, states, initial_state, nam, ndays, generator, warp_indicator,
                           base_size = 11, k = 1) {

  sol1 <- matrix(NA, nrow = ndays, ncol = m)
  count <- 0
  for (t in (1:ndays)/365) {
    count <- count + 1
    sol1[count, ] <- MJP_predict(m = m, s1 = c(initial_state), u = c(t),
      get_pars(nam), z = matrix(z[idx, ], nrow = 1),
      generator = generator, link_type_base = "exp", link_type_covs = "exp",
      covs_bin = TRUE, transient_dist_method = "eigen_decomp", warping = warp_indicator)
  }

  dpp <- as.data.frame(sol1)
  colnames(dpp) <- states
  dpp$time <- 1:ndays
  dpp_long <- tidyr::gather(dpp, key = "Column", value = "Probability", -time)
  dpp_long$Column <- factor(dpp_long$Column, levels = c("3","2B","2A","1","0"))

  ggplot(dpp_long, aes(x = time/365, y = Probability, color = Column)) +
  geom_line(size = 0.75 * k) +
  theme(
    text             = element_text(size = base_size * k, family = "serif"),
    axis.text.x      = element_text(size = rel(0.82), vjust = 0.3),
    axis.text.y      = element_text(size = rel(0.82)),
    legend.title     = element_blank(),
    legend.text      = element_text(size = rel(0.8)),
    legend.key.width = unit(0.8, "lines"),
    legend.key.height= unit(0.6, "lines"),
    panel.background = element_rect(fill = "white", color = "black"),
    panel.grid.minor = element_line(color = "lightgray"),
    legend.position  = c(0.55, 0.9),
    legend.direction = "horizontal",
    legend.background= element_rect(fill = "transparent", color = NA),
    legend.key       = element_rect(fill = "transparent", color = NA)
  ) +
  guides(color = guide_legend(nrow = 1)) +
  xlab("Time [years]") + ylab("Probability")
}

trans_prob_fig <- function(m, states, initial_state, nam, ndays, generator, warp_indicator,
                           base_size = 11, k = 1) {
  A1 <- matrix(0, m, m)
  for (i in 1:m) {
    A1[i, ] <- MJP_predict(m = m, s1 = c(i), u = c(ndays/365),
      get_pars(nam), z = matrix(z[idx, ], nrow = 1),
      generator = generator, link_type_base = "exp", link_type_covs = "exp",
      covs_bin = TRUE, transient_dist_method = "eigen_decomp", warping = warp_indicator)
  }
  colnames(A1) <- 1:5; rownames(A1) <- 1:5
  long <- reshape2::melt(A1)

  labs_map <- c("3","2B","2A","1","0")
  long$x <- factor(long$Var2, levels = 1:5, labels = labs_map)
  long$y <- factor(long$Var1, levels = 1:5, labels = labs_map)

  ggplot(long[long$value != 0, ], aes(x = x, y = y)) +
    geom_tile(aes(fill = value), linewidth = 0) +
    geom_text(aes(label = sprintf("%.5f", value)),
              color = "white", size = 3 * k, family = "serif") +
    scale_fill_gradient(low = "grey60", high = "black", guide = "none") +
    scale_x_discrete(drop = FALSE, expand = c(0,0)) +
    scale_y_discrete(drop = FALSE, expand = c(0,0), limits = rev(labs_map)) +
    coord_fixed() +
    labs(x = "To class", y = "From class") +
    theme(
      text             = element_text(size = base_size * k, family = "serif"),
      axis.text.x      = element_text(size = rel(0.82), vjust = 0.3),
      axis.text.y      = element_text(size = rel(0.82)),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.minor = element_blank()
    )
    

# helper: scale sizes based on device size versus a reference size (default 4x3 in)
size_scaler <- function(width_in, height_in, ref_w = 4, ref_h = 3) {
  # use the limiting dimension so proportions stay consistent
  min(width_in / ref_w, height_in / ref_h)
}

# trans_dist_fig = function(m, states, initial_state, nam, ndays, generator, warp_indicator){
	# sol1 = matrix(NA, nrow = ndays, ncol = m)
	# count = 0
	# for(t in (1:ndays)/365){
  		# count = count + 1
  		# sol1[count, ] = MJP_predict(m = m, s1 = c(initial_state), u = c(t), get_pars(nam), z = matrix(z[idx, ], nrow = 1), 		generator = generator, link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "eigen_decomp", warping = warp_indicator)
	# }
	# dpp <- as.data.frame(sol1)
	# colnames(dpp) <- states
	# dpp$time <- 1:ndays
	# dpp_long <- tidyr::gather(dpp, key = "Column", value = "Probability", -time)
	# dpp_long$Column <- factor(dpp_long$Column, levels = c("3", "2B", "2A", "1", "0")) 
	# p1 <- ggplot(data = dpp_long, aes(x = time/365, y = Probability, color = Column)) +
  	# geom_line(size = 0.75) +
  	# theme(
    # text = element_text(size = text.size, family = "serif"),
    # panel.background = element_rect(fill = "white", color = "black"),
    # panel.grid.minor = element_line(color = "lightgray"),
    # legend.position = c(0.55, 0.9),
    # legend.direction = "horizontal", # Set legend direction to horizontal
    # legend.background = element_rect(fill = "transparent", color = NA), # Set transparent background
    # legend.key = element_rect(fill = "transparent", color = NA) # Set transparent background for legend key
  	# ) + xlab("Time [years]") + ylab("Probability") + labs(color = "Classes")
  	# return(p1)
# }

# trans_prob_fig = function(m, states, initial_state, nam, ndays, generator, warp_indicator){
	# A1 = matrix(0, m,m)
	# for(i in 1:m){
		# A1[i,] = MJP_predict(m = m, s1 = c(i), u = c(ndays/365), get_pars(nam), z = matrix(z[idx, ], nrow = 1), generator = generator, link_type_base = "exp", link_type_covs = "exp", covs_bin = T, transient_dist_method = "eigen_decomp", warping = warp_indicator)
	# }
	# colnames(A1) = 1:5
	# rownames(A1) = 1:5
	# long <- melt(A1)  # Var1 = row (from), Var2 = col (to)
	

	# labs_map <- c("3","2B","2A","1","0")  # 1→"3", 2→"2B", 3→"2A", 4→"1", 5→"0"

	# # map indices to labels
	# long$x <- factor(long$Var2, levels = 1:5, labels = labs_map)  # columns = x
	# long$y <- factor(long$Var1, levels = 1:5, labels = labs_map)  # rows    = y

	# # plot: keep zeros as white space, but place cells exactly as matrix
	# p11 = ggplot(long[long$value != 0, ], aes(x = x, y = y)) +
  		# geom_tile(aes(fill = value), linewidth = 0) +
  		# geom_text(aes(label = sprintf("%.5f", value)),
            # color = "white", size = 3, family = "serif") +
  		# scale_fill_gradient(low = "grey60", high = "black", guide = "none") +
  		# scale_x_discrete(drop = FALSE, expand = c(0,0)) +
  		# scale_y_discrete(drop = FALSE, expand = c(0,0), limits = rev(labs_map)) +  # put row 1 on top
  		# coord_fixed() +
  		# labs(x = "To class", y = "From class") +
  		# theme(
    		# axis.text.x = element_text(size = 9, vjust = 0.3),
    		# axis.text.y = element_text(size = 9),
    		# text        = element_text(size = 11, family = "serif"),
    		# panel.background = element_rect(fill = "white", color = "black"),
    		# panel.grid.minor = element_blank()
  			# )
  	# return(p11)
# }

