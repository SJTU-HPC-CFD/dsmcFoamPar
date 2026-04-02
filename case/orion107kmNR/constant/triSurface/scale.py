import trimesh
import numpy as np
import math

def scale_stl_trimesh(input_file, output_file, scale_factor=1.5):
    """
    使用trimesh库读取和缩放STL文件，以中心点为缩放中心
    
    Args:
        input_file (str): 输入STL文件路径
        output_file (str): 输出STL文件路径
        scale_factor (float): 缩放因子，大于1放大，小于1缩小
    """
    # 加载STL文件
    mesh = trimesh.load_mesh(input_file)
    
    print(f"原始模型信息:")
    print(f"顶点数: {len(mesh.vertices)}")
    print(f"面数: {len(mesh.faces)}")
    print(f"边界框: {mesh.bounds}")
    print(f"中心点: {mesh.centroid}")
    
    # 获取中心点
    center = mesh.centroid
    print(f"\n缩放前中心点: {center}")
    
    # 创建缩放变换矩阵
    # 1. 先平移到原点
    translation_to_origin = trimesh.transformations.translation_matrix(-center)
    
    # 2. 创建缩放变换
    scale_matrix = trimesh.transformations.scale_matrix(
        scale_factor,  # 缩放因子
        [0, 0, 0]     # 缩放中心点（原点）
    )
    
    # 3. 平移回原位置
    translation_back = trimesh.transformations.translation_matrix(center)
    
    # 组合变换：平移到原点 -> 缩放 -> 平移回原位置
    combined_transform = np.dot(
        translation_back,
        np.dot(scale_matrix, translation_to_origin)
    )
    
    # 应用变换
    mesh.apply_transform(combined_transform)
    
    # 验证缩放后的中心点
    new_center = mesh.centroid
    print(f"缩放后中心点: {new_center}")
    print(f"中心点变化: {new_center - center}")
    
    # 计算缩放比例验证
    original_bounds = mesh.bounds
    scaled_bounds = mesh.bounds
    print(f"\n原始边界框尺寸: {original_bounds[1] - original_bounds[0]}")
    print(f"缩放后边界框尺寸: {scaled_bounds[1] - scaled_bounds[0]}")
    actual_scale = (scaled_bounds[1] - scaled_bounds[0]) / (original_bounds[1] - original_bounds[0])
    print(f"实际缩放比例: {actual_scale}")
    
    # 保存缩放后的模型
    mesh.export(output_file)
    print(f"\n缩放后的模型已保存到: {output_file}")
    
    return mesh

# 使用示例
if __name__ == "__main__":
    input_stl = "spline.stl"
    output_stl = "spline_new.stl"
    
    try:
        scaled_mesh = scale_stl_trimesh(
            input_stl,
            output_stl,
            scale_factor=1.5  # 放大1.5倍
        )
        
        print(f"\n最终模型信息:")
        print(f"顶点数: {len(scaled_mesh.vertices)}")
        print(f"面数: {len(scaled_mesh.faces)}")
        print(f"边界框: {scaled_mesh.bounds}")
        print(f"中心点: {scaled_mesh.centroid}")
        
    except Exception as e:
        print(f"错误: {e}")
