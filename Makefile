# Makefile for ahc_shc_anc_hao_with_rho project

# 编译器
FC = mpiifort

# 编译选项
FFLAGS = -CB -r8 -mkl

# 目标可执行文件
TARGET = ahc_shc_anc_hao_with_rho.x

# 源文件
SRC_MODULE = pauli.f90 fermi.f90
SRC_MAIN = ahc_shc_anc_hao_with_rho.f90

# 生成的对象文件
OBJ_MODULE = pauli.o fermi.o
OBJ_MAIN = ahc_shc_anc_hao_with_rho.o

# 默认目标
.PHONY: all
all: $(TARGET)

# 链接目标
$(TARGET): $(OBJ_MAIN) $(OBJ_MODULE)
	$(FC) $(FFLAGS) $(OBJ_MAIN) $(OBJ_MODULE) -o $(TARGET)

# 编译 pauli.f90 模块
$(OBJ_MODULE): $(SRC_MODULE)
	$(FC) $(FFLAGS) -c $(SRC_MODULE)

# 编译主程序
$(OBJ_MAIN): $(SRC_MAIN) $(OBJ_MODULE)
	$(FC) $(FFLAGS) -c $(SRC_MAIN)

# 清理规则
.PHONY: clean
clean:
	rm -f $(TARGET) $(OBJ_MODULE) $(OBJ_MAIN) pauli.mod
