# ===========================================================================
# Плоская задача упругости - квазистатическое нагружение 2D пластины с отверстием
# Наименование программы:   StressAnalysis
# Автор:                    D Klyukin
# Дата создания:            01.05.2026
# Версия:                   2.7
# Дата последней модификации: 07.06.2026
# ===========================================================================
# Описание:
#   Расчёт напряжений Мизеса в пластине с квадратным отверстием
#   при возрастающей вертикальной силе. Используется МКЭ на треугольной
#   сетке (RTriangle), решение СЛАУ с помощью LU-разложения,
#   визуализация в PNG с единой цветовой шкалой, сборка GIF.
# ===========================================================================
library(RTriangle)
library(Matrix)
library(fields)
library(gifski)

set.seed(42)

# --------------------- Параметры материала ---------------------------
E  <- 200e9 # Нержавеющая сталь
nu <- 0.3
D_mat <- E / (1 - nu^2) * matrix(c(1, nu, 0,
                                   nu, 1, 0,
                                   0, 0, (1 - nu) / 2), nrow = 3, byrow = TRUE)

# --------------------- Геометрия и сетка --------------------------
outer <- matrix(c(0,0, 1,0, 1,1, 0,1), ncol = 2, byrow = TRUE)
hole  <- matrix(c(0.4,0.4, 0.6,0.4, 0.6,0.6, 0.4,0.6), ncol = 2, byrow = TRUE)
all_vertices <- rbind(outer, hole)
edge_outer <- cbind(1:4, c(2:4, 1))
edge_hole  <- cbind(5:8, c(6:8, 5))   
all_edges  <- rbind(edge_outer, edge_hole)
hole_point <- matrix(c(0.5, 0.5), ncol = 2)

psl <- RTriangle::pslg(P = all_vertices, S = all_edges, H = hole_point)
tri <- RTriangle::triangulate(psl, a = 0.0001, q = 30)
nodes <- tri$P
tris  <- tri$T
N_nodes <- nrow(nodes)
N_tris  <- nrow(tris)

# --------------------- Сборка матрицы жёсткости -----------------
element_stiffness_2D <- function(i_tri) {
  idx <- tris[i_tri, ]
  xy <- nodes[idx, ]
  x1 <- xy[1,1]; y1 <- xy[1,2]
  x2 <- xy[2,1]; y2 <- xy[2,2]
  x3 <- xy[3,1]; y3 <- xy[3,2]
  area2 <- abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
  area <- area2 / 2.0
  if (area < 1e-12) return(matrix(0, 6, 6))
  b <- c(y2 - y3, y3 - y1, y1 - y2)
  c <- c(x3 - x2, x1 - x3, x2 - x1)
  B <- matrix(0, 3, 6)
  for (i in 1:3) {
    B[1, 2*i - 1] <- b[i]
    B[2, 2*i]     <- c[i]
    B[3, 2*i - 1] <- c[i]
    B[3, 2*i]     <- b[i]
  }
  B <- B / (2 * area)
  K_loc <- t(B) %*% D_mat %*% B * area
  return(K_loc)
}

# Глобальная матрица
K_i <- c(); K_j <- c(); K_val <- c()
for (t in 1:N_tris) {
  K_loc <- element_stiffness_2D(t)
  nd <- tris[t, ]
  dof <- as.vector(rbind(2*nd - 1, 2*nd))
  for (i in 1:6) {
    for (j in 1:6) {
      K_i   <- c(K_i, dof[i])
      K_j   <- c(K_j, dof[j])
      K_val <- c(K_val, K_loc[i, j])
    }
  }
}
K_global <- sparseMatrix(i = K_i, j = K_j, x = K_val, dims = c(2*N_nodes, 2*N_nodes))
K_global <- forceSymmetric(K_global)

# --------------------- Граничные условия (закрепления) ---------------
tol <- 1e-9
left_nodes   <- which(abs(nodes[,1]) < tol)
bottom_nodes <- which(abs(nodes[,2]) < tol)

for (node in left_nodes) {
  dof <- 2*node - 1
  K_global[dof, ] <- 0
  K_global[, dof] <- 0
  K_global[dof, dof] <- 1
}
for (node in bottom_nodes) {
  dof <- 2*node
  K_global[dof, ] <- 0
  K_global[, dof] <- 0
  K_global[dof, dof] <- 1
}

luK <- lu(K_global)

# --------------------- Параметры нагружения ---------------------------
force_levels <- seq(200, 2000, length.out = 10)
top_nodes <- which(abs(nodes[,2] - 1) < tol)
center_load <- which(nodes[top_nodes, 1] >= 0.3 & nodes[top_nodes, 1] <= 0.7)
load_nodes <- top_nodes[center_load]

# --------------------- ОПРЕДЕЛЕНИЕ МАКСИМАЛЬНОГО НАПРЯЖЕНИЯ ----------
cat("Определение максимального напряжения Мизеса при максимальной нагрузке...\n")
f_max <- max(force_levels)   # 2000 Н

F_vec_max <- rep(0, 2*N_nodes)
force_per_node_max <- f_max / length(load_nodes)
for (node in load_nodes) F_vec_max[2*node] <- F_vec_max[2*node] + force_per_node_max
F_vec_max[2*left_nodes - 1] <- 0
F_vec_max[2*bottom_nodes]   <- 0

U_max <- as.vector(solve(luK, F_vec_max))
u_max <- U_max[seq(1, 2*N_nodes, by = 2)]
v_max <- U_max[seq(2, 2*N_nodes, by = 2)]

sigma_vm_max <- 0
for (t in 1:N_tris) {
  idx <- tris[t, ]
  xy <- nodes[idx, ]
  x1 <- xy[1,1]; y1 <- xy[1,2]
  x2 <- xy[2,1]; y2 <- xy[2,2]
  x3 <- xy[3,1]; y3 <- xy[3,2]
  area2 <- abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
  area <- area2 / 2.0
  if (area < 1e-12) next
  b <- c(y2 - y3, y3 - y1, y1 - y2)
  c <- c(x3 - x2, x1 - x3, x2 - x1)
  B <- matrix(0, 3, 6)
  for (i in 1:3) {
    B[1, 2*i - 1] <- b[i]; B[2, 2*i] <- c[i]
    B[3, 2*i - 1] <- c[i]; B[3, 2*i] <- b[i]
  }
  B <- B / (2 * area)
  u_elem <- c(u_max[idx[1]], v_max[idx[1]], u_max[idx[2]], v_max[idx[2]], u_max[idx[3]], v_max[idx[3]])
  sigma <- D_mat %*% (B %*% u_elem)
  sigma_vm <- sqrt(sigma[1]^2 + sigma[2]^2 - sigma[1]*sigma[2] + 3*sigma[3]^2)
  if (sigma_vm > sigma_vm_max) sigma_vm_max <- sigma_vm
}
cat("Глобальный максимум напряжения Мизеса (при 2000 Н):", sigma_vm_max, "Па\n")

# -------------------- Подготовка палитры -------------------------
zlim <- c(0, sigma_vm_max)
pal <- colorRampPalette(c("blue", "cyan", "yellow", "red"))(100)
get_color <- function(val) {
  idx <- round((val - zlim[1]) / diff(zlim) * 99 + 1)
  idx <- pmax(1, pmin(100, idx))
  pal[idx]
}

# ---------------- Сохранение кадров для всех нагрузок -----------
graphics.off()

frames_dir <- "frames"
if (!dir.exists(frames_dir)) {
  dir.create(frames_dir, showWarnings = FALSE)
  if (!dir.exists(frames_dir)) {
    stop("Не удалось создать папку 'frames'. Проверьте права на запись.")
  }
}

abs_frames_dir <- file.path(getwd(), frames_dir)
cat("Кадры будут сохранены в:", abs_frames_dir, "\n")
cat("Сохранение кадров для всех уровней нагрузки...\n")

for (k in seq_along(force_levels)) {
  f_total <- force_levels[k]
  
  F_vec <- rep(0, 2*N_nodes)
  force_per_node <- f_total / length(load_nodes)
  for (node in load_nodes) F_vec[2*node] <- F_vec[2*node] + force_per_node
  F_vec[2*left_nodes - 1] <- 0
  F_vec[2*bottom_nodes]   <- 0
  
  U <- as.vector(solve(luK, F_vec))
  u <- U[seq(1, 2*N_nodes, by = 2)]
  v <- U[seq(2, 2*N_nodes, by = 2)]
  
  sigma_vm <- numeric(N_tris)
  for (t in 1:N_tris) {
    idx <- tris[t, ]
    xy <- nodes[idx, ]
    x1 <- xy[1,1]; y1 <- xy[1,2]
    x2 <- xy[2,1]; y2 <- xy[2,2]
    x3 <- xy[3,1]; y3 <- xy[3,2]
    area2 <- abs((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
    area <- area2 / 2.0
    if (area < 1e-12) next
    b <- c(y2 - y3, y3 - y1, y1 - y2)
    c <- c(x3 - x2, x1 - x3, x2 - x1)
    B <- matrix(0, 3, 6)
    for (i in 1:3) {
      B[1, 2*i - 1] <- b[i]; B[2, 2*i] <- c[i]
      B[3, 2*i - 1] <- c[i]; B[3, 2*i] <- b[i]
    }
    B <- B / (2 * area)
    u_elem <- c(u[idx[1]], v[idx[1]], u[idx[2]], v[idx[2]], u[idx[3]], v[idx[3]])
    sigma <- D_mat %*% (B %*% u_elem)
    sigma_vm[t] <- sqrt(sigma[1]^2 + sigma[2]^2 - sigma[1]*sigma[2] + 3*sigma[3]^2)
  }
  
  png_filename <- file.path(frames_dir, sprintf("frame_%04d.png", k))
  full_png_path <- file.path(getwd(), png_filename)
  
  success <- tryCatch({
    png(png_filename, width = 600, height = 500, bg = "white")
    plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "x", ylab = "y",
         main = paste0("Напряжение по Мизесу [Па] (сила = ", f_total, " Н)"))
    for (t in 1:N_tris) {
      polygon(nodes[tris[t,], 1], nodes[tris[t,], 2],
              col = get_color(sigma_vm[t]), border = NA)
    }
    lines(c(0,1,1,0,0), c(0,0,1,1,0), col = "red", lwd = 2)
    lines(c(0.4,0.6,0.6,0.4,0.4), c(0.4,0.4,0.6,0.6,0.4), col = "red", lwd = 2)
    fields::image.plot(zlim = zlim, col = pal, legend.only = TRUE,
                       smallplot = c(0.85,0.88, 0.2,0.8))
    dev.off()
    TRUE
  }, error = function(e) {
    message("Ошибка при сохранении кадра ", k, ": ", e$message)
    FALSE
  })
  
  if (success) {
    cat("Сохранён кадр", k, "/", length(force_levels), "->", full_png_path, "\n")
  } else {
    cat("Не удалось сохранить кадр", k, "\n")
  }
}

# --------------------- Создание GIF ----------------------------------
png_files <- list.files(frames_dir, pattern = "frame_.*\\.png", full.names = TRUE)
png_files <- sort(png_files)
cat("Найдено", length(png_files), "кадров в папке", abs_frames_dir, "\n")

if (length(png_files) == 0) {
  stop("Нет ни одного PNG-файла для создания GIF.")
}

gif_filename <- "stress_animation.gif"
gif_full_path <- file.path(getwd(), gif_filename)

gifski(png_files, 
       gif_file = gif_filename,
       width = 600,
       height = 500,
       delay = 0.5)

cat("GIF успешно создан:", gif_full_path, "\n")