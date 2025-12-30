#!/bin/bash
# ============================================
# MIRA 微服务部署脚本
# ============================================
#
# 使用方法:
#   ./deploy.sh [命令] [选项]
#
# 命令:
#   build     构建所有镜像
#   build-cpu 构建 CPU 版镜像
#   test      启动测试环境 (Gateway + MACE-ORB)
#   test-cpu  启动 CPU 测试环境 (无GPU)
#   up        启动所有服务
#   up-cpu    启动 CPU 模式所有服务
#   down      停止所有服务
#   logs      查看日志
#   status    查看服务状态
#   clean     清理未使用的镜像
#
# ============================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DOCKER_DIR="$PROJECT_ROOT/docker"

# 颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# 检查 Docker 环境
check_requirements() {
    log_info "检查环境..."
    
    if ! command -v docker &> /dev/null; then
        log_error "Docker 未安装"
        exit 1
    fi
    
    if ! command -v docker-compose &> /dev/null && ! docker compose version &> /dev/null; then
        log_error "Docker Compose 未安装"
        exit 1
    fi
}

# 检查 GPU 支持
check_gpu() {
    # 检查 NVIDIA Docker
    if docker run --rm --gpus all nvidia/cuda:12.1-base nvidia-smi &> /dev/null; then
        log_success "NVIDIA Docker 可用"
        nvidia-smi --query-gpu=index,name,memory.total --format=csv
        return 0
    else
        log_warn "NVIDIA Docker 不可用"
        return 1
    fi
}

# 构建 GPU 镜像
build_images() {
    log_info "构建 Docker 镜像 (GPU 版)..."
    cd "$PROJECT_ROOT"
    
    log_info "构建 Gateway..."
    docker build -f docker/Dockerfile.gateway -t mira-gateway:latest .
    
    log_info "构建 MACE-ORB Worker..."
    docker build -f docker/Dockerfile.mace-orb -t mira-mace-orb:latest .
    
    log_info "构建 FAIRChem-SevenNet Worker..."
    docker build -f docker/Dockerfile.fairchem-sevennet -t mira-fairchem-sevennet:latest .
    
    log_info "构建 MatGL Worker..."
    docker build -f docker/Dockerfile.matgl -t mira-matgl:latest .
    
    log_info "构建 GRACE Worker..."
    docker build -f docker/Dockerfile.grace -t mira-grace:latest .
    
    log_info "构建 MatterSim Worker..."
    docker build -f docker/Dockerfile.mattersim -t mira-mattersim:latest .
    
    log_success "所有 GPU 镜像构建完成"
    docker images | grep mira
}

# 构建 CPU 镜像
build_images_cpu() {
    log_info "构建 Docker 镜像 (CPU 版)..."
    cd "$PROJECT_ROOT"
    
    log_info "构建 Gateway..."
    docker build -f docker/Dockerfile.gateway -t mira-gateway:latest .
    
    log_info "构建 MACE-ORB Worker (CPU)..."
    docker build -f docker/Dockerfile.mace-orb-cpu -t mira-mace-orb-cpu:latest .
    
    log_info "构建 MatGL Worker (CPU)..."
    docker build -f docker/Dockerfile.matgl-cpu -t mira-matgl-cpu:latest .
    
    log_success "CPU 镜像构建完成"
    docker images | grep mira
}

# 启动 GPU 测试环境
start_test() {
    log_info "启动测试环境 (Gateway + MACE-ORB) - GPU..."
    cd "$DOCKER_DIR"
    
    docker compose -f docker-compose.test.yml up -d
    
    log_info "等待服务启动..."
    sleep 10
    
    log_info "检查服务状态..."
    docker compose -f docker-compose.test.yml ps
    
    log_success "测试环境已启动"
    log_info "访问 http://localhost:8000/docs 查看 API 文档"
}

# 启动 CPU 测试环境
start_test_cpu() {
    log_info "启动测试环境 (Gateway + MACE-ORB) - CPU..."
    cd "$DOCKER_DIR"
    
    docker compose -f docker-compose.cpu.yml up -d
    
    log_info "等待服务启动..."
    sleep 10
    
    log_info "检查服务状态..."
    docker compose -f docker-compose.cpu.yml ps
    
    log_success "CPU 测试环境已启动"
    log_info "访问 http://localhost:8000/docs 查看 API 文档"
    log_warn "注意: CPU 模式计算速度较慢，仅用于功能测试"
}

# 启动所有 GPU 服务
start_all() {
    log_info "启动所有微服务 (GPU)..."
    cd "$DOCKER_DIR"
    
    docker compose -f docker-compose.microservices.yml up -d
    
    log_info "等待服务启动..."
    sleep 30
    
    log_info "检查服务状态..."
    docker compose -f docker-compose.microservices.yml ps
    
    log_success "所有服务已启动"
    log_info "访问 http://localhost:8000/docs 查看 API 文档"
}

# 启动所有 CPU 服务
start_all_cpu() {
    log_info "启动所有微服务 (CPU 生产模式)..."
    cd "$DOCKER_DIR"
    
    docker compose -f docker-compose.cpu-prod.yml up -d
    
    log_info "等待服务启动 (CPU 模式需要较长时间)..."
    sleep 60
    
    log_info "检查服务状态..."
    docker compose -f docker-compose.cpu-prod.yml ps
    
    log_success "CPU 生产环境已启动"
    log_info "访问 http://localhost:8000/docs 查看 API 文档"
    log_warn "注意: CPU 模式计算速度较慢，适合小批量计算"
}

# 停止服务
stop_services() {
    log_info "停止服务..."
    cd "$DOCKER_DIR"
    
    docker compose -f docker-compose.microservices.yml down 2>/dev/null || true
    docker compose -f docker-compose.test.yml down 2>/dev/null || true
    docker compose -f docker-compose.cpu.yml down 2>/dev/null || true
    docker compose -f docker-compose.cpu-prod.yml down 2>/dev/null || true
    
    log_success "服务已停止"
}

# 查看日志
view_logs() {
    cd "$DOCKER_DIR"
    
    if [ -n "$2" ]; then
        docker compose -f docker-compose.microservices.yml logs -f "$2"
    else
        docker compose -f docker-compose.microservices.yml logs -f
    fi
}

# 查看状态
view_status() {
    log_info "服务状态:"
    cd "$DOCKER_DIR"
    
    docker compose -f docker-compose.microservices.yml ps 2>/dev/null || \
    docker compose -f docker-compose.test.yml ps 2>/dev/null || \
    log_warn "没有运行中的服务"
    
    echo ""
    log_info "健康检查:"
    curl -s http://localhost:8000/health 2>/dev/null | python3 -m json.tool || \
    log_warn "Gateway 未响应"
}

# 清理
clean_up() {
    log_info "清理未使用的 Docker 资源..."
    
    docker system prune -f
    docker image prune -f
    
    log_success "清理完成"
}

# 主函数
main() {
    case "$1" in
        build)
            check_requirements
            check_gpu
            build_images
            ;;
        build-cpu)
            check_requirements
            build_images_cpu
            ;;
        test)
            check_requirements
            check_gpu && start_test || start_test_cpu
            ;;
        test-cpu)
            check_requirements
            start_test_cpu
            ;;
        up|start)
            check_requirements
            check_gpu && start_all || start_all_cpu
            ;;
        up-cpu|start-cpu)
            check_requirements
            start_all_cpu
            ;;
        down|stop)
            stop_services
            ;;
        logs)
            view_logs "$@"
            ;;
        status|ps)
            view_status
            ;;
        clean)
            clean_up
            ;;
        *)
            echo "MIRA 微服务部署脚本"
            echo ""
            echo "使用方法: $0 [命令]"
            echo ""
            echo "GPU 命令:"
            echo "  build       构建所有 GPU Docker 镜像"
            echo "  test        启动 GPU 测试环境 (Gateway + MACE-ORB)"
            echo "  up          启动所有 GPU 服务"
            echo ""
            echo "CPU 命令 (无 GPU 环境使用):"
            echo "  build-cpu   构建 CPU Docker 镜像"
            echo "  test-cpu    启动 CPU 测试环境"
            echo "  up-cpu      启动 CPU 模式服务"
            echo ""
            echo "通用命令:"
            echo "  down        停止所有服务"
            echo "  logs        查看日志 (可选: logs [服务名])"
            echo "  status      查看服务状态"
            echo "  clean       清理未使用的 Docker 资源"
            echo ""
            echo "示例:"
            echo "  $0 build           # 构建 GPU 镜像"
            echo "  $0 build-cpu       # 构建 CPU 镜像 (无 GPU 服务器)"
            echo "  $0 test            # 自动检测: GPU 可用则用 GPU，否则用 CPU"
            echo "  $0 test-cpu        # 强制使用 CPU 模式"
            echo "  $0 up              # 启动生产环境 (自动检测)"
            echo "  $0 logs gateway    # 查看 Gateway 日志"
            ;;
    esac
}

main "$@"
