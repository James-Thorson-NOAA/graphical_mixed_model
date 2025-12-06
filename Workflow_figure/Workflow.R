
setwd( R'(C:\Users\James.Thorson\Desktop\Git\graphical_mixed_model\Workflow_figure)' )

library(igraph)

plot_graph <-
function( graph, label = names(V(graph)) ){           # rev(seq_len(ncol(A)))
  A = (as_adjacency_matrix(graph))   # [,rev(seq_len(nrow(A)))]
  image( as.matrix(A)[,rev(seq_len(ncol(A)))], x = seq_len(nrow(A)), y = seq_len(ncol(A)),
         col = c("white","black"), xaxt = "n", yaxt = "n",
         xlab = "", ylab = "" )
  axis(1, at = seq_len(nrow(A)), label = label )   # , las = 2
  axis(2, at = seq_len(nrow(A)), label = rev(label) )    # , las = 2
  box()
}


# Adjacency
G = make_graph( ~ s_1 - s_2 - s_3 - s_4 - s_5 - s_6 - s_7 - s_8 )
M = as_adjacency_matrix(G)
png( file = "adjacency.png", width = 4, height = 4, res = 200, units = "in" )
  par( mar = c(3,3,1,1), mgp = c(2,0.5,0), tck = -0.02 )
  #image(M, sub = "", xlab = "", ylab = "", col.regions = c("black") )
  plot_graph(t(G) )
dev.off()

#
G = make_empty_graph(3)
V(G)$name = c("c_1", "c_2", "c_3" )
G <- add_edges(G, c("c_3", "c_2"))
M = as_adjacency_matrix(G)
png( file = "predator_consumer.png", width = 4, height = 4, res = 200, units = "in" )
  par( mar = c(3,3,1,1), mgp = c(2,0.5,0), tck = -0.02 )
  plot_graph(t(G) )
dev.off()

# Adjacency
G = make_empty_graph(5)
V(G)$name = c("t_1", "t_2", "t_3", "t_4", "t_5" )
G <- add_edges(G, c("t_2", "t_1"))
G <- add_edges(G, c("t_3", "t_2"))
G <- add_edges(G, c("t_4", "t_3"))
G <- add_edges(G, c("t_5", "t_4"))
png( file = "lag.png", width = 4, height = 4, res = 200, units = "in" )
  par( mar = c(3,3,1,1), mgp = c(2,0.5,0), tck = -0.02 )
  #image(M, sub = "", xlab = "", ylab = "", col.regions = c("black") )
  plot_graph( t(G) )
dev.off()

# Time-identity
G = make_empty_graph(5)
V(G)$name = c("t_1", "t_2", "t_3", "t_4", "t_5" )
G <- add_edges(G, c("t_1", "t_1"))
G <- add_edges(G, c("t_2", "t_2"))
G <- add_edges(G, c("t_3", "t_3"))
G <- add_edges(G, c("t_4", "t_4"))
G <- add_edges(G, c("t_5", "t_5"))
png( file = "time_identity.png", width = 4, height = 4, res = 200, units = "in" )
  par( mar = c(3,3,1,1), mgp = c(2,0.5,0), tck = -0.02 )
  #image(M, sub = "", xlab = "", ylab = "", col.regions = c("black") )
  plot_graph( t(G) )
dev.off()
                                