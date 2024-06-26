#### CellTrek Functions ####
# traint ####
traint <- function (
    st_data, sc_data, st_assay='Spatial', sc_assay='scint', 
    norm='LogNormalize', nfeatures=2000, 
    cell_names='cell_names', coord_xy=c('imagerow', 'imagecol'), 
    gene_kept=NULL, ...) {
  st_data$id <- names(st_data$orig.ident)
  sc_data$id <- names(sc_data$orig.ident)
  sc_data$cell_names <- make.names(sc_data@meta.data[, cell_names])
  st_data$type <- 'st'
  sc_data$type <- 'sc'
  st_data$coord_x <- st_data@images[[1]]@coordinates[, coord_xy[1]]
  st_data$coord_y <- st_data@images[[1]]@coordinates[, coord_xy[2]]
  DefaultAssay(st_data) <- st_assay
  DefaultAssay(sc_data) <- sc_assay
  
  cat('Finding transfer anchors... \n')
  st_idx <- st_data$id
  sc_idx <- sc_data$id
  
  ## Integration features ##
  sc_st_list <- list(st_data=st_data, sc_data=sc_data)
  sc_st_features <- Seurat::SelectIntegrationFeatures(
    sc_st_list, nfeatures=nfeatures)
  if (!is.null(gene_kept)) {
    sc_st_features <- union(sc_st_features, gene_kept)
  }
  
  sc_st_features <- 
    sc_st_features[(sc_st_features %in% rownames(st_data[[st_assay]]@data)) & 
                     (sc_st_features %in% rownames(sc_data[[sc_assay]]@data))]
  cat('Using', length(sc_st_features), 'features for integration... \n')
  ###
  
  sc_st_anchors <- Seurat::FindTransferAnchors(
    reference = sc_data, query = st_data, 
    reference.assay = sc_assay, query.assay = st_assay,
    normalization.method = norm, features = sc_st_features, 
    reduction = 'cca')
  
  cat('Data transfering... \n')
  st_data_trans <- Seurat::TransferData(
    anchorset = sc_st_anchors, 
    refdata = GetAssayData(sc_data)[sc_st_features, ], 
    weight.reduction = 'cca')
  
  st_data@assays$transfer <- st_data_trans
  
  cat('Creating new Seurat object... \n')
  sc_st_meta <- dplyr::bind_rows(st_data@meta.data, sc_data@meta.data)
  counts_temp <- cbind(
    data.frame(st_data[['transfer']]@data),
    data.frame(sc_data[[sc_assay]]@data[sc_st_features, ] %>% data.frame))
  rownames(sc_st_meta) <- make.names(sc_st_meta$id)
  colnames(counts_temp) <- make.names(sc_st_meta$id)
  sc_st_int <- CreateSeuratObject(
    counts = counts_temp, 
    assay = 'traint', 
    meta.data = sc_st_meta, project = "SeuratProject"
  )
  
  cat('Scaling -> PCA -> UMAP... \n')
  sc_st_int <- ScaleData(
    sc_st_int, 
    features = sc_st_features, 
    assay = 'traint', layer = 'counts') %>%
    RunPCA(features = sc_st_features)
  sc_st_int <- RunUMAP(sc_st_int, dims = 1:30)
  sc_st_int@images <- st_data@images
  library(magrittr)
  sc_st_int@images[[1]]@coordinates <- data.frame(
    imagerow=sc_st_int@meta.data$coord_x,
    imagecol=sc_st_int@meta.data$coord_y) %>% 
    set_rownames(rownames(sc_st_int@meta.data))
  
  return (sc_st_int)
}

# celltrek ####
celltrek <- function (
    st_sc_int, int_assay='traint', sc_data=NULL, sc_assay='RNA', 
    reduction='pca', intp=T, intp_pnt=10000, intp_lin=F, nPCs=30, ntree=1000, 
    dist_thresh=.4, top_spot=10, spot_n=10, repel_r=5, repel_iter=10, 
    keep_model=F, ...) {
  
  dist_res <- celltrek_dist(
    st_sc_int=st_sc_int, int_assay=int_assay, reduction=reduction, intp=intp, 
    intp_pnt=intp_pnt, intp_lin=intp_lin, nPCs=nPCs, ntree=ntree, keep_model=T)
  
  spot_dis_intp <- mediaspot_dis_intp <- mediaspot_dis_intp <- median(
    unlist(dbscan::kNN(dist_res$coord_df[, c('coord_x', 'coord_y')], k=4)$dist))
  
  if (is.null(repel_r)) {repel_r=spot_dis_intp/4}
  
  sc_coord_list <- celltrek_chart(
    dist_mat=dist_res$celltrek_dist, coord_df=dist_res$coord_df, 
    dist_cut=ntree*dist_thresh, top_spot=top_spot, spot_n=spot_n, 
    repel_r=repel_r, repel_iter=repel_iter)
  
  sc_coord_raw <- sc_coord_list[[1]]
  sc_coord <- sc_coord_list[[2]]
  
  cat('Creating Seurat Object... \n')
  if (!is.null(sc_data)) {
    cat('sc data...')
    sc_data$id <- Seurat::Cells(sc_data)
    
    sc_out <- CreateSeuratObject(
      counts=sc_data[[sc_assay]]@data[, sc_coord$id_raw] %>% 
        set_colnames(sc_coord$id_new), 
      project='celltrek', 
      assay=sc_assay, 
      meta.data=sc_data@meta.data[sc_coord$id_raw, ] %>% 
        dplyr::rename(id_raw=id) %>% 
        mutate(id_new=sc_coord$id_new) %>% 
        set_rownames(sc_coord$id_new))
    
    sc_out@meta.data <- dplyr::left_join(
      sc_out@meta.data, sc_coord) %>% 
      data.frame %>% 
      set_rownames(sc_out$id_new)
    
    sc_coord_raw_df <- CreateDimReducObject(
      embeddings=sc_coord_raw %>% 
        dplyr::mutate(
          coord1=coord_y, coord2=max(coord_x)+min(coord_x)-coord_x
        ) %>%
        dplyr::select(c(coord1, coord2)) %>% 
        set_rownames(sc_coord_raw$id_new) %>% 
        as.matrix, assay=sc_assay, key='celltrek_raw')
    
    sc_coord_dr <- CreateDimReducObject(
      embeddings=sc_coord %>% 
        dplyr::mutate(
          coord1=coord_y, coord2=max(coord_x)+min(coord_x)-coord_x) %>% 
        dplyr::select(c(coord1, coord2)) %>% 
        set_rownames(sc_coord$id_new) %>% 
        as.matrix, assay=sc_assay, key='celltrek')
    
    sc_out@reductions$celltrek <- sc_coord_dr
    sc_out@reductions$celltrek_raw <- sc_coord_raw_df
    if ('pca' %in% names(sc_data@reductions)) {
      sc_pca_dr <- CreateDimReducObject(
        embeddings=sc_data@reductions$pca@cell.embeddings[sc_coord$id_raw, ] %>% 
          set_rownames(sc_coord$id_new) %>% as.matrix, assay=sc_assay, key='pca'
      )
      sc_out@reductions$pca <- sc_pca_dr
    }
    if ('umap' %in% names(sc_data@reductions)) {
      sc_umap_dr <- CreateDimReducObject(
        embeddings=sc_data@reductions$umap@cell.embeddings[sc_coord$id_raw,] %>% 
          set_rownames(sc_coord$id_new) %>% as.matrix, assay=sc_assay, 
        key='umap'
      )
      sc_out@reductions$umap <- sc_umap_dr
    }
    if ('tsne' %in% names(sc_data@reductions)) {
      sc_tsne_dr <- CreateDimReducObject(
        embeddings=sc_data@reductions$tsne@cell.embeddings[sc_coord$id_raw, ] %>% 
          set_rownames(sc_coord$id_new) %>% 
          as.matrix, assay=sc_assay, key='tsne')
      sc_out@reductions$tsne <- sc_tsne_dr
    }
  } else {
    cat('no sc data...')
    sc_out <- CreateSeuratObject(
      counts=st_sc_int[[int_assay]]@data[, sc_coord$id_raw] %>% 
        set_colnames(sc_coord$id_new), 
      project='celltrek', assay=int_assay, 
      meta.data=st_sc_int@meta.data[sc_coord$id_raw, ] %>%
        dplyr::rename(id_raw=id) %>%
        mutate(id_new=sc_coord$id_new) %>%
        set_rownames(sc_coord$id_new)
    )
    sc_out$coord_x <- sc_coord$coord_x[match(sc_coord$id_new, sc_out$id_new)]
    sc_out$coord_y <- sc_coord$coord_y[match(sc_coord$id_new, sc_out$id_new)]
    
    sc_out[[int_assay]]@scale.data <- 
      st_sc_int[[int_assay]]@scale.data[, sc_coord$id_raw] %>% 
      set_colnames(sc_coord$id_new)
    sc_coord_raw_df <- CreateDimReducObject(
      embeddings=sc_coord_raw %>% 
        dplyr::mutate(
          coord1=coord_y, coord2=max(coord_x)+min(coord_x)-coord_x) %>%
        dplyr::select(c(coord1, coord2)) %>% 
        set_rownames(sc_coord_raw$id_new) %>% 
        as.matrix, assay=sc_assay, key='celltrek_raw'
    )
    sc_coord_dr <- CreateDimReducObject(
      embeddings = sc_coord %>% 
        dplyr::mutate(
          coord1=coord_y, coord2=max(coord_x)+min(coord_x)-coord_x) %>% 
        dplyr::select(c(coord1, coord2)) %>% 
        set_rownames(sc_coord$id_new) %>% 
        as.matrix, assay=int_assay, key='celltrek'
    )
    sc_out@reductions$celltrek <- sc_coord_dr
    sc_out@reductions$celltrek_raw <- sc_coord_raw_df
    if ('pca' %in% names(st_sc_int@reductions)) {
      sc_pca_dr <- CreateDimReducObject(
        embeddings=st_sc_int@reductions$pca@cell.embeddings[sc_coord$id_raw, ] %>% 
          set_rownames(sc_coord$id_new) %>%
          as.matrix, assay=int_assay, key='pca'
      )
      sc_out@reductions$pca <- sc_pca_dr
    }
    if ('umap' %in% names(st_sc_int@reductions)) {
      sc_umap_dr <- CreateDimReducObject(
        embeddings=st_sc_int@reductions$umap@cell.embeddings[sc_coord$id_raw,] %>% 
          set_rownames(sc_coord$id_new) %>%
          as.matrix, assay=int_assay, key='umap'
      )
      sc_out@reductions$umap <- sc_umap_dr
    }
    if ('tsne' %in% names(st_sc_int@reductions)) {
      sc_tsne_dr <- CreateDimReducObject(
        embeddings=st_sc_int@reductions$tsne@cell.embeddings[sc_coord$id_raw,] %>%
          set_rownames(sc_coord$id_new) %>%
          as.matrix, assay=int_assay, key='tsne'
      )
      sc_out@reductions$tsne <- sc_tsne_dr
    }
  }
  sc_out@images <- st_sc_int@images
  sc_out@images[[1]]@assay <- DefaultAssay(sc_out)
  sc_out@images[[1]]@coordinates <- data.frame(
    imagerow=sc_coord$coord_x, imagecol=sc_coord$coord_y) %>%
    set_rownames(sc_coord$id_new)
  sc_out@images[[1]]@scale.factors$spot_dis <- dist_res$spot_d
  sc_out@images[[1]]@scale.factors$spot_dis_intp <- spot_dis_intp
  
  output <- list(celltrek=sc_out)
  if (keep_model) {
    output[[length(output)+1]] <- dist_res$model
    names(output)[length(output)] <- 'model'
  }
  return(output)
}

# celltrek_dist ####
celltrek_dist <- function (
    st_sc_int, int_assay='traint', reduction='pca', intp = T, intp_pnt=10000, 
    intp_lin=F, nPCs=30, ntree=1000, keep_model=T) {
  DefaultAssay(st_sc_int) <- int_assay
  kNN_dist <- dbscan::kNN(
    na.omit(st_sc_int@meta.data[, c('coord_x', 'coord_y')]), k=6)$dist
  spot_dis <- median(kNN_dist) %>% round
  cat('Distance between spots is:', spot_dis, '\n')
  
  st_sc_int$id <- names(st_sc_int$orig.ident)
  st_idx <- st_sc_int$id[st_sc_int$type=='st']
  sc_idx <- st_sc_int$id[st_sc_int$type=='sc']
  meta_df <- data.frame(st_sc_int@meta.data)
  
  st_sc_int_pca <- 
    st_sc_int@reductions[[reduction]]@cell.embeddings[, 1:nPCs] %>% 
    data.frame %>%
    mutate(id=st_sc_int$id, type=st_sc_int$type, class=st_sc_int$cell_names,
           coord_x=st_sc_int$coord_x, coord_y=st_sc_int$coord_y)
  st_pca <- st_sc_int_pca %>% 
    dplyr::filter(type=='st') %>% 
    dplyr::select(-c(id:class))
  
  ## Interpolation ##
  ## Uniform sampling ##
  if (intp) {
    cat ('Interpolating...\n')
    spot_ratio <- intp_pnt/nrow(st_pca)
    st_intp_df <- apply(st_pca[, c('coord_x', 'coord_y')], 1, function(row_x) {
      runif_test <- runif(1)
      if (runif_test < spot_ratio%%1) {
        theta <- runif(ceiling(spot_ratio), 0, 2*pi)
        alpha <- sqrt(runif(ceiling(spot_ratio), 0, 1))
        coord_x <- row_x[1] + (spot_dis/2)*sin(theta)*alpha
        coord_y <- row_x[2] + (spot_dis/2)*cos(theta)*alpha
      } else {
        theta <- runif(floor(spot_ratio), 0, 2*pi)
        alpha <- sqrt(runif(floor(spot_ratio), 0, 1))
        coord_x <- row_x[1] + (spot_dis/2)*sin(theta)*alpha
        coord_y <- row_x[2] + (spot_dis/2)*cos(theta)*alpha
      }
      data.frame(coord_x, coord_y)
    }) %>% Reduce(rbind, .)
    
    st_intp_df <- apply(st_pca[, 1:nPCs], 2, function(col_x) {
      akima::interpp(x=st_pca$coord_x, y=st_pca$coord_y, z=col_x,
                     linear=intp_lin, xo=st_intp_df$coord_x, 
                     yo=st_intp_df$coord_y) %>%
        magrittr::extract2('z')
    }) %>% data.frame(., id='X', type='st_intp', st_intp_df) %>% na.omit
    st_intp_df$id <- make.names(st_intp_df$id, unique = T)
    st_sc_int_pca <- bind_rows(st_sc_int_pca, st_intp_df)
  }
  
  cat('Random Forest training... \n')
  ## Training on ST ##
  data_train <- st_sc_int_pca %>% 
    dplyr::filter(type=='st') %>% 
    dplyr::select(-c(id:class))
  rf_train <- randomForestSRC::rfsrc(
    Multivar(coord_x, coord_y) ~ ., data_train, block.size=5, ntree=ntree)
  
  cat('Random Forest prediction...  \n')
  ## Testing on all ##
  data_test <- st_sc_int_pca
  rf_pred <- randomForestSRC::predict.rfsrc(
    rf_train, newdata=data_test[, c(1:nPCs)], distance='all')
  
  cat('Making distance matrix... \n')
  rf_pred_dist <- 
    rf_pred$distance[data_test$type=='sc', data_test$type!='sc'] %>%
    set_rownames(data_test$id[data_test$type=='sc']) %>% 
    set_colnames(data_test$id[data_test$type!='sc'])
  
  output <- list()
  output$spot_d <- spot_dis
  output$celltrek_dist <- rf_pred_dist
  output$coord_df <- st_sc_int_pca[, c('id', 'type', 'coord_x', 'coord_y')] %>%
    dplyr::filter(type!='sc') %>% 
    magrittr::set_rownames(.$id) %>% 
    dplyr::select(-id)
  if (keep_model) {
    output$model <- rf_train
  }
  return (output)
}

# celltrek_chart ####
celltrek_chart <- function (dist_mat, coord_df, dist_cut=500, top_spot=10, 
                            spot_n=10, repel_r=5, repel_iter=10) {
  cat('Making graph... \n')
  dist_mat[dist_mat>dist_cut] <- NA
  dist_mat_dt <- data.table::data.table(Var1=rownames(dist_mat), dist_mat)
  dist_edge_list <- data.table::melt(dist_mat_dt, id=1, na.rm=T)
  colnames(dist_edge_list) <- c('Var1', 'Var2', 'value')
  dist_edge_list$val_rsc <- scales::rescale(
    dist_edge_list$value, to=c(0, repel_r))
  dist_edge_list$Var1 %<>% as.character
  dist_edge_list$Var2 %<>% as.character
  dist_edge_list$Var1_type <- 'sc'
  dist_edge_list$Var2_type <- 'non-sc'
  
  cat('Pruning graph...\n')
  dist_edge_list_sub <- dplyr::inner_join(
    dist_edge_list %>% 
      group_by(Var1) %>% 
      top_n(n=top_spot, wt=-value), 
    dist_edge_list %>% group_by(Var2) %>% 
      top_n(n=spot_n, wt=-value)) %>% 
    data.frame
  
  cat('Spatial Charting SC data...\n')
  sc_coord <- sc_coord_raw <- data.frame(
    id_raw=dist_edge_list_sub$Var1, id_new=make.names(dist_edge_list_sub$Var1, 
                                                      unique = T))
  sc_coord$coord_x <- 
    sc_coord_raw$coord_x <- 
    coord_df$coord_x[match(dist_edge_list_sub$Var2, rownames(coord_df))]
  sc_coord$coord_y <- 
    sc_coord_raw$coord_y <- 
    coord_df$coord_y[match(dist_edge_list_sub$Var2, rownames(coord_df))]
  ## Add noise ##
  theta <- runif(nrow(dist_edge_list_sub), 0, 2*pi)
  alpha <- sqrt(runif(nrow(dist_edge_list_sub), 0, 1))
  sc_coord$coord_x <- 
    sc_coord$coord_x + dist_edge_list_sub$val_rsc*sin(theta)*alpha
  sc_coord$coord_y <- 
    sc_coord$coord_y + dist_edge_list_sub$val_rsc*cos(theta)*alpha
  ## Point repelling ##
  cat('Repelling points...\n')
  sc_repel_input <- data.frame(
    sc_coord[, c('coord_x', 'coord_y')], repel_r=repel_r)
  sc_repel <- packcircles::circleRepelLayout(
    sc_repel_input, sizetype='radius', maxiter=repel_iter)
  sc_coord$coord_x <- sc_repel$layout$x
  sc_coord$coord_y <- sc_repel$layout$y
  return(list(sc_coord_raw, sc_coord))
}

# celltrek_vis ####
celltrek_vis <- function (celltrek_df, img, scale_fac) 
{
  img_fact <- scale_fac
  img_temp <- img
  img_data <- celltrek_df
  img_data$coord_x_new = img_data$coord_y * img_fact
  img_data$coord_y_new = dim(img_temp)[1] - img_data$coord_x * 
    img_fact
  if (!("id_new" %in% colnames(img_data))) 
    img_data$id_new <- rownames(img_data)
  app <- list(ui = fluidPage(titlePanel("CellTrek visualization"), 
                             sidebarLayout(sidebarPanel(width = 3, selectInput(inputId = "color_inp", 
                                                                               label = "Color", choices = colnames(img_data), selected = "None"), 
                                                        radioButtons(inputId = "color_typ", label = "Type", 
                                                                     choices = c("Categorical", "Continuous"), selected = "Categorical"), 
                                                        selectInput(inputId = "shape_inp", label = "Shape", 
                                                                    choices = c("None", colnames(img_data)), selected = "None"), 
                                                        actionButton("Plot", "Plot"), tags$hr(), textInput("colID", 
                                                                                                           "Add Type:"), actionButton("AddID", "Add"), 
                                                        downloadButton("downloadData", "Download"), tags$hr(), 
                                                        actionButton("StopID", "Stop")), mainPanel(plotlyOutput("CellTrek", 
                                                                                                                height = "1000px", width = "1200px"), DT::DTOutput("Tab_temp")))), 
              server = function(input, output, session) {
                options(warn = -1)
                data_react <- reactiveValues()
                data_react$DF <- data.frame(id_new = character(), 
                                            coord_x = numeric(), coord_y = numeric(), add_col = character())
                observeEvent(input$StopID, {
                  stopApp()
                })
                observeEvent(input$Plot, {
                  color_var <- isolate(input$color_inp)
                  type_var <- isolate(input$color_typ)
                  shape_var <- isolate(input$shape_inp)
                  if (type_var == "Categorical") {
                    output$CellTrek <- renderPlotly({
                      img_data$color_var <- factor(img_data[, 
                                                            color_var])
                      if (shape_var == "None") {
                        img_data$shape_var <- ""
                      } else {
                        img_data$shape_var <- factor(img_data[, 
                                                              shape_var])
                      }
                      if (length(levels(img_data$color_var)) <= 
                          9) {
                        pnt_colors <- brewer.pal(length(levels(img_data$color_var)), 
                                                 "Set1")
                      } else {
                        pnt_colors <- colorRampPalette(brewer.pal(9, 
                                                                  "Set1"))(length(levels(img_data$color_var)))
                      }
                      plotly::plot_ly(d = img_data, x = ~coord_x_new, 
                                      y = ~coord_y_new, customdata = ~id_new, 
                                      color = ~color_var, type = "scatter", 
                                      mode = "markers", text = ~color_var, symbol = ~shape_var, 
                                      colors = pnt_colors, marker = list(line = list(color = "rgb(1, 1, 1)", 
                                                                                     width = 0.5), size = 8, opacity = 0.8)) %>% 
                        plotly::layout(xaxis = list(range = c(0, 
                                                              dim(img_temp)[2]), showgrid = FALSE, 
                                                    showline = FALSE), yaxis = list(range = c(0, 
                                                                                              dim(img_temp)[1]), showgrid = FALSE, 
                                                                                    showline = FALSE), images = list(source = plotly::raster2uri(as.raster(img_temp)), 
                                                                                                                     x = 0, y = 0, sizex = dim(img_temp)[2], 
                                                                                                                     sizey = dim(img_temp)[1], xref = "x", 
                                                                                                                     yref = "y", xanchor = "left", yanchor = "bottom", 
                                                                                                                     layer = "below", sizing = "stretch"))
                    })
                  }
                  if (type_var == "Continuous") {
                    output$CellTrek <- renderPlotly({
                      img_data$color_var <- img_data[, color_var]
                      if (shape_var == "None") {
                        img_data$shape_var <- ""
                      } else {
                        img_data$shape_var <- factor(img_data[, 
                                                              shape_var])
                      }
                      plotly::plot_ly(d = img_data, x = ~coord_x_new, 
                                      y = ~coord_y_new, customdata = ~id_new, 
                                      color = ~color_var, type = "scatter", 
                                      mode = "markers", text = ~color_var, symbol = ~shape_var, 
                                      colors = c("#377EB8", "white", "#E41A1C"), 
                                      marker = list(line = list(color = "rgb(1, 1, 1)", 
                                                                width = 0.5), size = 8, opacity = 0.8)) %>% 
                        plotly::layout(xaxis = list(range = c(0, 
                                                              dim(img_temp)[2]), showgrid = FALSE, 
                                                    showline = FALSE), yaxis = list(range = c(0, 
                                                                                              dim(img_temp)[1]), showgrid = FALSE, 
                                                                                    showline = FALSE), images = list(source = plotly::raster2uri(as.raster(img_temp)), 
                                                                                                                     x = 0, y = 0, sizex = dim(img_temp)[2], 
                                                                                                                     sizey = dim(img_temp)[1], xref = "x", 
                                                                                                                     yref = "y", xanchor = "left", yanchor = "bottom", 
                                                                                                                     layer = "below", sizing = "stretch"))
                    })
                  }
                })
                observeEvent(input$AddID, {
                  tab_temp <- event_data("plotly_selected")
                  if (!is.null(tab_temp)) {
                    data_temp <- img_data[match(tab_temp$customdata, 
                                                img_data$id_new), c("id_new", "coord_x", 
                                                                    "coord_y")]
                    data_temp$add_col <- input$colID
                    data_react$DF <- bind_rows(data_react$DF, 
                                               data_temp)
                    output$Tab_temp <- renderDataTable(data_react$DF)
                  }
                })
                output$downloadData <- downloadHandler(filename = function() {
                  paste("changemyname.csv", sep = ",")
                }, content = function(file) {
                  write.csv(data_react$DF, file, row.names = F)
                })
              })
  shiny::runApp(app)
}
