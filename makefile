# コンパイラ設定
CC = icx
FC = ifx  # Fortranコンパイラ
CFLAGS = -O0 -Wall -Wextra -Isrc/c -Isrc/c/jacobs
FFLAGS = -O0 -fpe3 -check all -warn all # Fortran用のオプション
# CFLAGS = -g -O0 -Isrc/c -Isrc/c/jacobs
# FFLAGS = -g -O0 -fpe3
 # Fortran用のオプション
LDFLAGS = -qmkl  # 数学ライブラリをリンク

# # コンパイラ設定
# CC = gcc
# FC = gfortran  # Fortranコンパイラ
# CFLAGS = -O2 -Wall -Wextra -Isrc/c -Isrc/c/jacobs
# FFLAGS = -O2 -Wall -Wextra  # Fortran用のオプション
# LDFLAGS = -llapack -lblas -lm  # 数学ライブラリをリンク

# ソースファイルの場所
C_SRC_DIR = src/c
JACOBS_SRC_DIR = src/c/jacobs
F_SRC_DIR = src/fortran
OBJ_DIR = obj
BIN_DIR = bin

# Cソースファイル一覧
SRCS_C = $(wildcard $(C_SRC_DIR)/*.c)

# Cソースファイル一覧
SRCS_JACOBS = $(wildcard $(JACOBS_SRC_DIR)/*.c)

# Fortranソースファイル一覧
SRCS_F90 = $(wildcard $(F_SRC_DIR)/*.f90)

# ヘッダファイル一覧（オプションで管理）
DEPS = $(wildcard include/*.h)

# Cのオブジェクトファイル
OBJS_C = $(patsubst $(C_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS_C))

# Cのオブジェクトファイル
OBJS_JACOBS = $(patsubst $(JACOBS_SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS_JACOBS))

# Fortranのオブジェクトファイル
OBJS_F90 = $(patsubst $(F_SRC_DIR)/%.f90, $(OBJ_DIR)/%.o, $(SRCS_F90))

# 実行ファイルの名前
TARGET = $(BIN_DIR)/main

# ディレクトリ作成
$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

# すべてのオブジェクトファイルをリンクして実行ファイルを作成
$(TARGET): $(OBJS_C) $(OBJS_JACOBS) $(OBJS_F90)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# 各Cファイルをコンパイル
$(OBJ_DIR)/%.o: $(C_SRC_DIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

# 各Cファイルをコンパイル
$(OBJ_DIR)/%.o: $(JACOBS_SRC_DIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

# 各Fortranファイルをコンパイル
$(OBJ_DIR)/%.o: $(F_SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c -o $@ $<

$(OBJ_DIR)/main.o: $(OBJ_DIR)/pv_parallel.o $(OBJ_DIR)/cema.o
$(OBJ_DIR)/cema.o: $(OBJ_DIR)/pv_parallel.o

# クリーンアップ
.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) *.mod src/fortran/*.mod *__genmod.f90
