CC = gcc
CFLAGS = -O2 -Wall -Wextra -Iinclude
LDFLAGS = -lm  # 数学ライブラリをリンク

# ソースファイルの場所
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# ソースファイル一覧
SRCS = $(wildcard $(SRC_DIR)/*.c)

# ヘッダファイル一覧（オプションで管理）
DEPS = $(wildcard include/*.h)

# オブジェクトファイル
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

# 実行ファイルの名前
TARGET = $(BIN_DIR)/main

# ディレクトリ作成
$(shell mkdir -p $(OBJ_DIR) $(BIN_DIR))

# すべてのオブジェクトファイルをリンクして実行ファイルを作成
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# 各Cファイルをコンパイル
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

# クリーンアップ
.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
