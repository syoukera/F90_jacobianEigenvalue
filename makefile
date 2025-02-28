# コンパイラ設定
CC = gcc
FC = gfortran  # Fortranコンパイラ
CFLAGS = -O2 -Wall -Wextra
FFLAGS = -O2 -Wall -Wextra  # Fortran用のオプション
LDFLAGS = -llapack -lblas -lm  # 数学ライブラリをリンク

# ソースファイルの場所
C_SRC_DIR = src/c
F_SRC_DIR = src/fortran
OBJ_DIR = obj
BIN_DIR = bin

# Cソースファイル一覧
SRCS_C = $(wildcard $(C_SRC_DIR)/*.c)

# Fortranソースファイル一覧
SRCS_F90 = $(wildcard $(F_SRC_DIR)/*.f90)

# ヘッダファイル一覧（オプションで管理）
DEPS = $(wildcard include/*.h)

# Cのオブジェクトファイル
OBJS_C = $(patsubst $(C_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS_C))

# Fortranのオブジェクトファイル
OBJS_F90 = $(patsubst $(F_SRC_DIR)/%.f90, $(OBJ_DIR)/%.o, $(SRCS_F90))

# 実行ファイルの名前
TARGET = $(BIN_DIR)/main

# ディレクトリ作成
$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

# すべてのオブジェクトファイルをリンクして実行ファイルを作成
$(TARGET): $(OBJS_C) $(OBJS_F90)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# 各Cファイルをコンパイル
$(OBJ_DIR)/%.o: $(C_SRC_DIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

# 各Fortranファイルをコンパイル
$(OBJ_DIR)/%.o: $(F_SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJ_DIR)/main.o: $(OBJ_DIR)/pv_parallel.o $(OBJ_DIR)/cema.o

# クリーンアップ
.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) *.mod src/fortran/*.mod
