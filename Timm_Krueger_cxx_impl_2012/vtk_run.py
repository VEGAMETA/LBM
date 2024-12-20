import os, re
import pyvista as pv
import imageio

# Папки с данными
fluid_dir = "vtk_fluid"
particle_dir = "vtk_particle"
def sort_key(filename):
    match = re.search(r"t(\d+)", filename)
    return int(match.group(1)) if match else 0
# Список файлов
fluid_files = sorted([os.path.join(fluid_dir, f) for f in os.listdir(fluid_dir) if f.endswith(".vtk")], key=sort_key)
particle_files = sorted([os.path.join(particle_dir, f) for f in os.listdir(particle_dir) if f.endswith(".vtk")], key=sort_key)

# Проверка на одинаковое количество файлов
assert len(fluid_files) == len(particle_files), "Количество файлов в папках должно совпадать"

# Настройки визуализации
plotter = pv.Plotter(off_screen=True)
frames = []

# Генерация кадров
for fluid_file, particle_file in zip(fluid_files, particle_files):
    # Загрузка данных
    fluid_data = pv.RectilinearGrid(fluid_file)
    particle_data = pv.PolyData(particle_file)

    # Очистка старой сцены
    plotter.clear()

    # Добавление данных
    plotter.add_mesh(fluid_data, scalars="density_difference", cmap="viridis", opacity=0.8)
    plotter.add_mesh(particle_data, color="red", point_size=5)

    # Настройка камеры (опционально)
    plotter.camera_position = [(50, 20, 135), (50, 20, 0), (0, -1, 0)]

    # Сохранение кадра
    img = plotter.screenshot(return_img=True)
    frames.append(img)

# Создание видео
output_video = "simulation.mp4"
imageio.mimsave(output_video, frames, fps=10)
print(f"Видео сохранено: {output_video}")
