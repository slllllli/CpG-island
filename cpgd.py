import tkinter as tk
from tkinter import ttk, messagebox, simpledialog
import math
import random
import threading
#markov model
CPGMATRIX = {
    'A': {'A': 0.180, 'C': 0.274, 'G': 0.426, 'T': 0.120},
    'C': {'A': 0.171, 'C': 0.368, 'G': 0.274, 'T': 0.188},
    'G': {'A': 0.161, 'C': 0.339, 'G': 0.375, 'T': 0.125},
    'T': {'A': 0.079, 'C': 0.355, 'G': 0.384, 'T': 0.182}
}

NOCPGMATRIX = {
    'A': {'A': 0.300, 'C': 0.205, 'G': 0.285, 'T': 0.210},
    'C': {'A': 0.322, 'C': 0.298, 'G': 0.078, 'T': 0.302},
    'G': {'A': 0.248, 'C': 0.246, 'G': 0.298, 'T': 0.208},
    'T': {'A': 0.177, 'C': 0.239, 'G': 0.292, 'T': 0.292}
}

def analyze_cpg(sequence):
    lg_cpg=0
    lg_no_cpg=0
    if len(sequence)<2:
        messagebox.showerror("错误", "序列长度必须大于等于2")
    for i in range(len(sequence)-1):
        pre_nuc=sequence[i].upper()
        curr_nuc=sequence[i+1].upper()
        lg_cpg+=math.log(CPGMATRIX[pre_nuc][curr_nuc])
        lg_no_cpg+=math.log(NOCPGMATRIX[pre_nuc][curr_nuc])
    lg=lg_cpg-lg_no_cpg

    return lg, "分析完成"

class CpGApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("CpG岛判断器 (马尔可夫链模型)")
        self.geometry('600x600')

        self.style = ttk.Style(self)
        self.style.theme_use('default')
        self.threshold = 0

        self.create_widgets()

    def create_widgets(self):
        #顶部菜单栏
        menubar = tk.Menu(self)
        self.config(menu=menubar)
        model_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="模型参数", menu=model_menu)
        model_menu.add_command(label="转移概率矩阵", command=self.show_matrix_window)
        model_menu.add_command(label="阈值设定", command=self.set_threshold)

        
        #主框架
        main_frame = ttk.Frame(self, padding="5")#主框架
        main_frame.pack(fill=tk.BOTH, expand=True)
        #状态栏
        self.status_label = ttk.Label(main_frame, text="就绪", relief=tk.SUNKEN, anchor=tk.W)
        self.status_label.pack(fill=tk.X, side=tk.TOP, pady=(5, 5))
        # 在状态栏后添加阈值显示
        threshold_frame = ttk.Frame(main_frame)
        threshold_frame.pack(fill=tk.X, pady=(5, 0))
        self.threshold_label = ttk.Label(threshold_frame, text=f"当前阈值: {self.threshold}（阈值为对数似然值分界点）")
        self.threshold_label.pack(side=tk.LEFT)

        # --- 输入区域 ---
        input_frame = ttk.LabelFrame(main_frame, text="DNA序列输入", padding="10")#输入区域、内边距
        input_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 5))#填充、扩展、外边距
        # --- 输入框 ---
        self.input_text = tk.Text(input_frame, height=5, font=("Consolas", 12))#输入框，初始高度
        self.input_text.pack(fill=tk.BOTH, expand=True)
        # --- 按钮区域 ---
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.BOTH)

        ttk.Button(button_frame, text="随机生成序列", command=self.generate_random_sequence).pack(fill=tk.X,pady=(5, 5))
        ttk.Button(button_frame, text="开始分析", command=self.start_analysis_thread).pack(fill=tk.X)

    def start_analysis_thread(self):
        sequence = self.input_text.get("1.0", tk.END).strip()
        if not sequence:
            messagebox.showwarning("输入为空", "请输入DNA序列后再进行判断。")
            return

        valid_chars = {'A', 'C', 'G', 'T'}
        if not all(char.upper() in valid_chars for char in sequence):
            messagebox.showerror("输入错误", "序列包含无效字符！请只输入 A, C, G, T。")
            return

        self.status_label.config(text="正在分析中...")

        # 在新线程中执行耗时操作，避免界面卡顿
        thread = threading.Thread(target=self.run_analysis, args=(sequence,))
        thread.daemon = True
        thread.start()

    def run_analysis(self, sequence): # 在新线程中执行耗时操作
        log_odds, status = analyze_cpg(sequence)
        # 分析完成后，在主线程中更新UI
        self.after(0, self.show_result_dialog, log_odds)

    def show_result_dialog(self, log_odds):
        self.status_label.config(text="分析完成")
        # 使用Sigmoid函数将log_odds转换为(0, 1)之间的概率
        # S(x) = 1 / (1 + e^(-x))
        lg_new = log_odds-self.threshold
        try:
            cpg_confidence_raw = 1 / (1 + math.exp(-lg_new))
        except OverflowError:
            cpg_confidence_raw = 1.0 if lg_new > 0 else 0.0

        cpg_confidence_percent = cpg_confidence_raw * 100
        non_cpg_confidence_percent = (1 - cpg_confidence_raw) * 100

        # --- 创建结果弹窗 ---
        result_window = tk.Toplevel(self)
        result_window.title("分析结果")
        result_window.geometry("450x300")
        result_window.transient(self)
        result_window.grab_set()

        # 居中显示
        result_window.update_idletasks()
        x = (result_window.winfo_screenwidth() // 2) - (result_window.winfo_width() // 2)
        y = (result_window.winfo_screenheight() // 2) - (result_window.winfo_height() // 2)
        result_window.geometry(f"+{x}+{y}")

        # --- 结果内容 ---
        main_frame = ttk.Frame(result_window, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)

        # 对数似然比
        ttk.Label(main_frame, text="对数似然比:", font=("Helvetica", 12)).grid(row=0, column=0, sticky=tk.W, pady=5)
        ttk.Label(main_frame, text=f"{log_odds:.4f}", font=("Helvetica", 12, "bold")).grid(row=0, column=1, sticky=tk.W, pady=5)

        # 最终判断
        final_judgment = "是 CpG 岛" if lg_new > 0 else "非 CpG 岛"
        judgment_color = "green" if lg_new > 0 else "red"
        ttk.Label(main_frame, text="最终判断:", font=("Helvetica", 12)).grid(row=1, column=0, sticky=tk.W, pady=5)
        ttk.Label(main_frame, text=final_judgment, font=("Helvetica", 14, "bold"), foreground=judgment_color).grid(row=1, column=1, sticky=tk.W, pady=5)

        # --- 置信度显示区域 ---
        ttk.Separator(main_frame, orient='horizontal').grid(row=2, column=0, columnspan=2, sticky="ew", pady=15)
        ttk.Label(main_frame, text="模型置信度:", font=("Helvetica", 12, "bold")).grid(row=3, column=0, columnspan=2, sticky=tk.W, pady=5)

        # CpG岛置信度
        ttk.Label(main_frame, text="是 CpG 岛:", font=("Helvetica", 11)).grid(row=4, column=0, sticky=tk.W, padx=(20, 5), pady=2)
        cpg_progressbar = ttk.Progressbar(main_frame, length=200, mode='determinate')
        cpg_progressbar.grid(row=4, column=1, sticky=tk.W, pady=2)
        cpg_progressbar['value'] = cpg_confidence_percent
        cpg_label = ttk.Label(main_frame, text=f"{cpg_confidence_percent:.2f}%", font=("Helvetica", 11))
        cpg_label.grid(row=4, column=2, sticky=tk.W, padx=(5, 0), pady=2)

        # 非CpG岛置信度
        ttk.Label(main_frame, text="非 CpG 岛:", font=("Helvetica", 11)).grid(row=5, column=0, sticky=tk.W, padx=(20, 5), pady=2)
        non_cpg_progressbar = ttk.Progressbar(main_frame, length=200, mode='determinate')
        non_cpg_progressbar.grid(row=5, column=1, sticky=tk.W, pady=2)
        non_cpg_progressbar['value'] = non_cpg_confidence_percent
        non_cpg_label = ttk.Label(main_frame, text=f"{non_cpg_confidence_percent:.2f}%", font=("Helvetica", 11))
        non_cpg_label.grid(row=5, column=2, sticky=tk.W, padx=(5, 0), pady=2)
        
        # --- 关闭按钮 ---
        ttk.Button(main_frame, text="关闭", command=result_window.destroy).grid(row=6, column=0, columnspan=3, pady=20)

    def generate_random_sequence(self):
        length_str = simpledialog.askstring("随机序列", "请输入要生成的序列长度 (例如: 200):", parent=self)
        if length_str and length_str.isdigit():
            length = int(length_str)
            if length > 0:
                random_sequence = "".join(random.choices('ACGT', k=length))
                self.input_text.delete("1.0", tk.END)
                self.input_text.insert("1.0", random_sequence)
                self.status_label.config(text=f"已生成 {length} bp 的随机序列。")
            else:
                messagebox.showerror("错误", "长度必须大于0。")

    def show_matrix_window(self):
        matrix_window = tk.Toplevel(self)
        matrix_window.title("模型参数：转移概率矩阵")
        matrix_window.geometry("650x400")
        matrix_window.transient(self)
        
        # 居中显示
        matrix_window.update_idletasks()
        x = (matrix_window.winfo_screenwidth() // 2) - (matrix_window.winfo_width() // 2)
        y = (matrix_window.winfo_screenheight() // 2) - (matrix_window.winfo_height() // 2)
        matrix_window.geometry(f"+{x}+{y}")

        main_frame = ttk.Frame(matrix_window, padding="10")
        main_frame.pack(fill=tk.BOTH, expand=True)

        main_frame.rowconfigure(0, weight=1)
        main_frame.rowconfigure(1, weight=1)
        main_frame.columnconfigure(0, weight=1)

        cpg_frame = ttk.LabelFrame(main_frame, text="CpG 岛模型转移概率", padding="10")
        cpg_frame.grid(row=0, column=0, sticky="nsew") # padx 在右侧留出一点间距
        self.create_matrix_table(cpg_frame, CPGMATRIX)

        non_cpg_frame = ttk.LabelFrame(main_frame, text="非 CpG 岛模型转移概率", padding="10")
        non_cpg_frame.grid(row=1, column=0, sticky="nsew") # padx 在左侧留出一点间距
        self.create_matrix_table(non_cpg_frame, NOCPGMATRIX)
    def set_threshold(self):
        threshold_str = simpledialog.askstring("阈值设置", f"当前阈值: {self.threshold}\n请输入新的阈值 (例如: 0.5):", parent=self)
        if threshold_str and threshold_str.replace('.', '', 1).isdigit():
            threshold = float(threshold_str)
        self.threshold = threshold
        self.threshold_label.config(text=f"当前阈值: {self.threshold}")
        messagebox.showinfo("成功", f"阈值已更新为: {self.threshold}")

    def create_matrix_table(self, parent, matrix_data):
        # 创建Treeview作为表格
        headers = ["Prev \\ Curr"] + list(matrix_data['A'].keys())
        tree = ttk.Treeview(parent, columns=headers, show='headings')
        
        # 设置列标题
        for col in headers:
            tree.heading(col, text=col)
            tree.column(col, anchor=tk.CENTER, width=100)

        # 添加数据行
        for prev_char, row_data in matrix_data.items():
            # 格式化数据为字符串，保留3位小数
            formatted_values = [f"{val:.3f}" for val in row_data.values()]
            tree.insert('', tk.END, values=[prev_char] + formatted_values)
        
        # 添加滚动条
        scrollbar = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=tree.yview)
        tree.configure(yscrollcommand=scrollbar.set)
        
        tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)


# --- 4. 程序入口 ---
if __name__ == "__main__":
    app = CpGApp()
    app.mainloop()