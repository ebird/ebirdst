library(ebirdst)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(terra)
library(rnaturalearth)

# download data for Asian House-Martin (Delichon dasypus)
# イワツバメ (Delichon dasypus) のデータをダウンロード
ebirdst_download_status("Delichon dasypus")

# load weekly raster data
# 週次ラスターデータをロードする
abd_weekly <- load_raster("ashmar1")

# load seasonal raster data
# 季節的なラスターデータをロードする
abd_seasonal <- load_raster("ashmar1", period = "seasonal")

# spatial boundaries for Akita and Fukuoka prefectures
# 秋田県と福岡県の空間的境界
boundaries <- ne_states(country = "Japan") |>
  filter(name %in% c("Akita", "Fukuoka")) |>
  select(prefecture = name) |>
  st_transform(crs = crs(abd_weekly))
plot(boundaries)

# extract values within prefecture boundaries and calculate average
# 都道府県境界内の値を抽出し、平均を算出する
chronology <- extract(abd_weekly, boundaries,
                      fun = "mean", na.rm = TRUE,
                      bind = TRUE, ID = FALSE)
# convert to a data frame
# データフレームに変換する
chronology <- as.data.frame(chronology)

# transform to wide format for plotting
# プロット用にワイドフォーマットに変換する
chronology <- pivot_longer(
  chronology,
  cols = starts_with("X"),
  names_to = "week", values_to = "abd",
  names_transform = \(x) as.Date(x, format = "X%Y.%m.%d")
)

# 移住の年代順をプロットする
ggplot(chronology) +
  aes(x = week, y = abd, color = prefecture) +
  geom_line() +
  geom_point() +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(x = "Week",
       y = "Mean relative abundance",
       title = "Migration chronology for Delichon dasypus")
