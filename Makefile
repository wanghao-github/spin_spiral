# Makefile for ahc_shc_anc_hao_with_rho project

# 编译器
FC = mpiifort

# 编译选项
FFLAGS = -CB -r8 -mkl

# 目标可执行文件
TARGET = ahc_shc_anc_hao_with_rho.x
SPIRAL_TARGET = spiral.x

# 源文件
SRC_MODULE = pauli.f90 fermi.f90
SRC_MAIN = ahc_shc_anc_hao_with_rho.f90
SRC_SPIRAL = anc_hao21.f90

# 生成的对象文件
OBJ_MODULE = pauli.o fermi.o
OBJ_MAIN = ahc_shc_anc_hao_with_rho.o
OBJ_SPIRAL = anc_hao21.o

# 生成的模块文件
MOD_FILES = pauli.mod fermi.mod

# 默认目标
.PHONY: all
all: $(TARGET) $(SPIRAL_TARGET)

# 链接 ahc_shc_anc_hao_with_rho 目标
$(TARGET): $(OBJ_MAIN) $(OBJ_MODULE)
	$(FC) $(FFLAGS) $(OBJ_MAIN) $(OBJ_MODULE) -o $(TARGET)

# 编译模块文件（依赖关系：先编译生成 .mod 文件）
$(OBJ_MODULE): %.o: %.f90
	$(FC) $(FFLAGS) -c $<

# 编译 ahc_shc_anc_hao_with_rho 主程序文件（依赖模块）
$(OBJ_MAIN): $(SRC_MAIN) $(MOD_FILES)
	$(FC) $(FFLAGS) -c $(SRC_MAIN)

# 编译 anc_hao21 模块文件（依赖模块）
$(OBJ_SPIRAL): $(SRC_SPIRAL) $(MOD_FILES)
	$(FC) $(FFLAGS) -c $(SRC_SPIRAL) -o $(OBJ_SPIRAL)

# 链接 spiral 目标
$(SPIRAL_TARGET): $(OBJ_SPIRAL) $(OBJ_MODULE)
	$(FC) $(FFLAGS) $(OBJ_SPIRAL) $(OBJ_MODULE) -o $(SPIRAL_TARGET)

# 清理规则
.PHONY: clean
clean:
	rm -f $(TARGET) $(SPIRAL_TARGET) $(OBJ_MODULE) $(OBJ_MAIN) $(OBJ_SPIRAL) *.mod
