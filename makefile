# コンパイラ設定
CC = gcc
FC = gfortran  # Fortranコンパイラ
CFLAGS = -O2 -Wall -Wextra -Iinclude
FFLAGS = -O2 -Wall -Wextra  # Fortran用のオプション
LDFLAGS = -llapack -lblas -lm  # 数学ライブラリをリンク

# ソースファイルの場所
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Cソースファイル一覧
SRCS_C = $(wildcard $(SRC_DIR)/*.c)

# Fortranソースファイル一覧
SRCS_F90 = $(wildcard $(SRC_DIR)/*.f90)

# ヘッダファイル一覧（オプションで管理）
DEPS = $(wildcard include/*.h)

# Cのオブジェクトファイル
OBJS_C = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS_C))

# Fortranのオブジェクトファイル
OBJS_F90 = $(patsubst $(SRC_DIR)/%.f90, $(OBJ_DIR)/%.o, $(SRCS_F90))

# 実行ファイルの名前
TARGET = $(BIN_DIR)/main

# ディレクトリ作成
$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

# すべてのオブジェクトファイルをリンクして実行ファイルを作成
$(TARGET): $(OBJS_C) $(OBJS_F90)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# 各Cファイルをコンパイル
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

# 各Fortranファイルをコンパイル
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

# クリーンアップ
.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
