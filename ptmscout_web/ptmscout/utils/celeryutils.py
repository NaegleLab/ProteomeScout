from celery.result import GroupResult


def get_percent_complete(import_process_id):
    res = GroupResult(import_process_id)
    children = float(len(res.children))
    completed = float(res.completed_count())
    
    return completed / children