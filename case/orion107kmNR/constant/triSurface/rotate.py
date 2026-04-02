import trimesh
import numpy as np
import math

def rotate_stl_trimesh(input_file, output_file, angle_degrees=45):
    """
    使用trimesh库读取和旋转STL文件
    """
    # 加载STL文件
    mesh = trimesh.load_mesh(input_file)
    
    print(f"原始模型信息:")
    print(f"顶点数: {len(mesh.vertices)}")
    print(f"面数: {len(mesh.faces)}")
    print(f"边界框: {mesh.bounds}")
    print(f"中心点: {mesh.centroid}")
    
    # 创建绕Y轴的旋转变换
    angle_rad = math.radians(angle_degrees)
    rotation_transform = trimesh.transformations.rotation_matrix(
        angle_rad, [0, 1, 0]  # Y轴
    )
    
    # 应用旋转变换
    mesh.apply_transform(rotation_transform)
    
    # 计算当前中心点
    current_center = mesh.centroid
    print(f"旋转后中心点: {current_center}")
    
    # 创建平移变换，将中心点移动到z=1.5
    # 计算需要的平移向量：目标z=1.5，保持x和y坐标不变
    target_z = 1.5
    translation_vector = [0, 0, target_z - current_center[2]]
    translation_transform = trimesh.transformations.translation_matrix(translation_vector)
    
    # 应用平移变换
    mesh.apply_transform(translation_transform)
    
    # 保存旋转和平移后的模型
    mesh.export(output_file)
    print(f"旋转并平移后的模型已保存到: {output_file}")
    print(f"最终中心点: {mesh.centroid}")
    
    return mesh

# 使用示例
if __name__ == "__main__":
    input_stl = "orion.stl"
    output_stl = "orion_new.stl"
    
    try:
        rotated_mesh = rotate_stl_trimesh(input_stl, output_stl, 45)
        
        print(f"\n最终模型信息:")
        print(f"顶点数: {len(rotated_mesh.vertices)}")
        print(f"面数: {len(rotated_mesh.faces)}")
        print(f"边界框: {rotated_mesh.bounds}")
        print(f"中心点: {rotated_mesh.centroid}")
        
    except Exception as e:
        print(f"错误: {e}")